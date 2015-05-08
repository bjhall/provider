#!/usr/bin/perl -w
use strict;
use threads;
use Getopt::Long;
use Data::Dumper;


my $SAMTOOLS_PATH = "samtools";
my $DEFAULT_SNPBED = "HPA_1000G_final_38.bed";
my $DEFAULT_XYBED  = "xy_38.bed";

my %opt = &get_options;
my %meta_data = &get_bampaths( $opt{bam}, $opt{nocheck} );
my ($variant_data, $chromosomes)  = &read_bed( $opt{bed} );

my ($xy_data, $xy_chromosomes);
if ( $opt{ bedxy } ) {
    ($xy_data, $xy_chromosomes) = &read_bed( $opt{bedxy} );
    foreach my $chr (keys %{$xy_chromosomes}) {
	$chromosomes->{$chr}++;
    }
}

unless ( $opt{ nocheck } or &matching_chr_names( (keys %meta_data)[0], $opt{bed} ) ) {
    error( "Chromosome names do not match between BED and BAMs", 1 );
}

my %data = &get_base_freqs_from_bams( \%meta_data, $opt{bed}, $opt{bedxy}, $chromosomes );

my %sample = &determine_sex( \%data, $xy_data );
&do_genotyping( \%data, $variant_data );

&print_genotype_table( \%data, \%sample, $opt{out}, \%meta_data, $variant_data );

&detect_unexpected( \%data, \%meta_data, \%sample );





########################

sub detect_unexpected{
    my( $data, $meta_data, $sample_data ) = @_;

    my (%dist, %seen);
    foreach my $samp1 ( keys %$meta_data ) {
	my $name1 = $meta_data->{$samp1}->{name};
	$seen{$samp1} = 1;

	my $anno_sex = $meta_data->{$samp1}->{sex};
	my $pred_sex = $sample_data->{$name1}->{sex}->{'Y:2844149'};
	print "UNEXPECTED SEX: $name1 (Predicted:$pred_sex, Annotated:$anno_sex)\n" if ($anno_sex ne "-" and $anno_sex ne $pred_sex);
	

	foreach my $samp2 (keys %$meta_data) {
	    next if $seen{$samp2};

	    my $name2 = $meta_data->{$samp2}->{name};
	    unless (defined( $dist{$name2}->{$name1} )) {
		my $dist = distance( $data, $name1, $name2 );
		$dist{$name1}->{$name2} = $dist;
		my ($ind1, $ind2) = ( $meta_data->{$samp1}->{individual}, $meta_data->{$samp2}->{individual} );

		if ($ind1 eq $ind2 and $dist > 0.05) {
		    printf "UNEXPECTED DIFFERENT: %s - %s (%.2f%%)\n", $name1, $name2, 100*(1-$dist{$name1}->{$name2});
		}
		elsif ($ind1 ne $ind2 and $dist < 0.05) {
		    printf "UNEXPECTED IDENTICAL: %s - %s (%.2f%%)\n", $name1, $name2, 100*(1-$dist{$name1}->{$name2});
		}
		
	    }
	}
    }
}

sub distance{
    my( $data, $id1, $id2 ) = @_;

    my ( $identical, $tot ) = ( 0, 0 );
    foreach my $loc (keys %$data) {
	if( defined($data->{$loc}->{samples}->{$id1}->{gt}) and defined($data->{$loc}->{samples}->{$id2}->{gt}) ) {
	    $identical++ if $data->{$loc}->{samples}->{$id1}->{gt} eq $data->{$loc}->{samples}->{$id2}->{gt};
	    $tot++;
	}
    }

    return ($tot-$identical) / $tot;
}

sub read_bed{
    my $bed = shift;
    open( BED, $bed );
    my (%var, %chr);
    while( <BED> ) {
	chomp;
	my @a = split /\t/;
	$chr{$a[0]} ++;
	if ($a[3]) {
	    $var{$a[0].":".$a[2]} = $a[3];
	}
	else {
	    $var{$a[0].":".$a[2]} = $a[0].":".$a[2];
	}
    }
    return( \%var, \%chr );
}

sub print_genotype_table{
    my( $snp_data, $sample_data, $out, $annotation, $variants ) = @_;

    open( GT, ">".$out.".genotypes" );

    if( ! $opt{'long'} ) {
	print GT "loc";
	print GT "\t$_" foreach sort keys %{ ((values %$snp_data)[0])->{samples} };
	print GT "\n";
    }
  
    my ($tot_callable, $tot_sites, $HW_pass, $tot_loc);
    foreach my $loc ( sort keys %$snp_data ) {
	my $snp_id = $opt{position} ? $loc : $variants->{$loc};

	print GT $snp_id if ! $opt{long};

	foreach my $sid ( sort keys %{ $snp_data->{$loc}->{samples} } ) {
	    my $genotype = $snp_data->{$loc}->{samples}->{$sid}->{gt};

	    if( $opt{long} ) {
		print GT "$sid\t$sid\t".$snp_id."\t". 
		         plink_gt( $snp_data->{$loc}->{samples}->{$sid}->{basecall} ) ."\n";
	    }
	    else {
		print GT "\t". ( defined( $genotype ) ? $genotype : "NA" ) if ! $opt{long};
	    }
	}

	print GT "\n" if ! $opt{long};
	$tot_callable += $snp_data->{$loc}->{Ncallable};
	$tot_sites    += $snp_data->{$loc}->{Ntotal};
	$HW_pass      += $snp_data->{$loc}->{HW};
	$tot_loc      ++;
    }   
    close GT;

    open (STATS, ">".$out.".stats");
    printf STATS "Average callability: %.2f%%\nSites in H-W equilibrium: %d / %d\n", 100*($tot_callable/$tot_sites),
                                                                                     $HW_pass, $tot_loc;
    close STATS;

    # Output file with sex prediction
    open (SAMPLES, ">".$out.".sex");
    print SAMPLES "sample\tsex\n";
    foreach my $bam ( sort keys %$annotation ) {
	my $sid = $annotation->{$bam}->{name};
	print SAMPLES $sid;
	foreach my $loc ( sort keys %$xy_data ) {
	    print SAMPLES "\t". $sample_data->{$sid}->{sex}->{$loc}."\t".($annotation->{$bam}->{sex} or "-");
	}
	print SAMPLES "\n";
    }
    close SAMPLES;
}


sub plink_gt{
    my $bc = shift;
    return "0\t0" if !defined($bc) or $bc eq "unclear" or $bc eq "lowdata" ;
    return "$bc\t$bc" if length($bc) == 1;
    my ($a1, $a2) = split /\//, $bc;
    return "$a1\t$a2";
}
 

sub determine_sex{
    my( $data, $loci ) = @_;

    my %sample;
    foreach my $loc ( keys %$loci ) {
	foreach my $sid (sort keys %{ $data->{$loc}->{samples} }) {
	    my $depth = $data->{$loc}->{samples}->{$sid}->{depth};
	    
	    # FIXME: Arbitrary numbers...
	    if ($depth >= 50) {
		$sample{$sid}->{sex}->{$loc} = "M";
	    }
	    elsif ($depth <= 5) {
		$sample{$sid}->{sex}->{$loc} = "F";
	    }
	    else {
		$sample{$sid}->{sex}->{$loc} = "unclear";
	    }

	}
	delete( $data->{$loc} );
    }
    return %sample;
}


sub do_genotyping{
    my( $data, $loci ) = @_;

    foreach my $loc ( keys %$loci ) {
	my $num_samples;

	# SNP loci: Do base calling
	my %gt_cnt;
	foreach my $sid (sort keys %{ $data->{$loc}->{samples} }) {
    
	    my %f = %{ $data->{$loc}->{samples}->{$sid}->{bases} };
	    my $tot = $data->{$loc}->{samples}->{$sid}->{depth};
	    my ($max, $second) = ( large(\%f, 1), large(\%f, 2) );

	    # Require at least 6 reads in the SNP position.
	    my $gt;
	    if ($f{$max} >= 6) {
		my $ap = ($f{$max} - $f{$second}) / $tot;
	
		if ($ap < 0.6) { # Heterozygote
		    $gt = join "/", sort( $max,$second );
		}
		elsif ($ap > 0.9) { # Homozygote
		    $gt = $max;
		}
		else {
		    $gt = 'unclear';
		}
	    }
	    else {
		$gt = 'lowdata';
	    }
	    $data{$loc}->{samples}->{$sid}->{basecall} = $gt;
	    $gt_cnt{$gt}++ if $gt ne "unclear" and $gt ne "lowdata";
	    $num_samples++;
	}
	$data->{$loc}->{Ntotal} = $num_samples;

	# Skipping tri-allelic sites
	if (keys %gt_cnt > 3) {
	    print STDERR "Skipping $loc: More than more than 2 alleles detected!\n";
	    next;
	}

	# Assign genotypes to 0,1,2 (0 = common homozygote, 1 = heterozygote, 2 = rare homozygote)
	my $first = 1;
	my %num;
	foreach my $gtype (sort {$gt_cnt{$b} <=> $gt_cnt{$a}} keys %gt_cnt) {
	    if ($gtype =~ /\//) {
		$num{$gtype} = 1;
	    }
	    else {
		$num{$gtype} = 0 if $first;
		$num{$gtype} = 2 if !$first;
		$first = 0;
	    }
	    $gt_cnt{ $num{$gtype} } = $gt_cnt{ $gtype };
	}

	# Add 0,1,2 genotypes to hash
	foreach my $sid (keys %{$data->{$loc}->{samples}}) {
	    my $call = $data->{$loc}->{samples}->{$sid}->{basecall};
	    if ($call ne "unclear" and $call ne "lowdata") {
		$data->{$loc}->{samples}->{$sid}->{gt} = $num{ $call }
	    }
	}

	# Check for Hardy-Weinberg equilibrium
	my ($AA, $AB, $BB) = ( ($gt_cnt{0} or 0), ($gt_cnt{1} or 0), ($gt_cnt{2} or 0) );
	$data{$loc}->{HW}        = HW_chi2test( $AA, $AB, $BB );
	$data{$loc}->{Ncallable} = $AA + $AB + $BB;
	$data{$loc}->{MAF}       = ( $AB + 2 * $BB ) / 2 * ($AA+$AB+$BB);
    }   
}


# Run mpileup for a single chromosome in a thread.
sub mpileup_in_thread {
    my ($chr, $out, $bams) = @_;
    my $id = threads->tid();
    system( $SAMTOOLS_PATH." mpileup -l ". $opt{out}.".tmp.bed -r $chr $bams".
	    " > $out.$chr 2> $opt{out}.mpileup.$chr.log" );
    threads->exit();
}


# Run mpilup on bams and count number of A, C, G and Ts in each position.
sub get_base_freqs_from_bams{
    my( $file_data, $snp_fn, $xy_fn, $chromosomes ) = @_;

    my $pile_out = $opt{out}.".pile.tmp";

    merge_files( [$snp_fn, $xy_fn], $opt{out}.".tmp.bed" );

    my @threads;
    if( $opt{overwrite} or !-s $pile_out ) {

	# Run mpileup in threads, split by chromosomes
	if( $opt{thread} ) {
	    my $i = 0;
	    my @pileup_files;
	    foreach my $chr ( sort { $chromosomes->{$b} <=> $chromosomes->{$a} } keys %{$chromosomes} ) {
		$threads[$i] = threads->create( \&mpileup_in_thread, ( $chr, $pile_out, join( " ", sort keys %$file_data) ) );
		push (@pileup_files, "$pile_out.$chr");
		$i++;
	    }
	    $_->join() for @threads;

	    merge_files( \@pileup_files, $pile_out );
	    unlink $_ foreach @pileup_files;
	}

	# Run mpileup unthreaded
	else {
	    system( $SAMTOOLS_PATH." mpileup -l ". $opt{out}.".tmp.bed ". join( " ", sort keys %$file_data ).
		    " > $pile_out 2> $opt{out}.mpileup.log" );
	}

	print STDERR " Done\n";
    }
    else {
	print STDERR "WARNING: Prior pileup file found, resuming with that file. Use --overwrite to override\n";
    }

    my %frq;

    open( PILE, $pile_out );
    while (<PILE>) {
	my @p = split /\t/;
	my $loc = $p[0].":".$p[1];
	
	my $i;
	foreach my $id ( sort keys %$file_data ) {
	    $i += 3;
	    my ($depth, $base, $qual) = ( $p[$i], $p[$i+1], $p[$i+2] );
	    my $bases = count_bases( $base, $qual );
	    $frq{$loc}->{samples}->{ $file_data->{$id}->{name} }->{bases} = $bases;
	    $frq{$loc}->{samples}->{ $file_data->{$id}->{name} }->{depth} = sum( values %{ $bases } );
	}
    }

    return %frq
}


sub merge_files1{ # FIXME: REMOVE
    my( $a, $b, $out ) = @_;

    my( @a, @b );

    open( A, $a );
    @a = <A>;
    close A;

    if ($b) {
	open( B, $b );
	@b = <B>;
	close B;
    }

    open( OUT, ">$out" );
    print OUT join( "", @a );
    print OUT join( "", @b ) if $b;
    close OUT;
}

sub merge_files{
    my ($files, $out) = @_;
    
    open( OUT, ">$out" );
    foreach (@$files) {
	open (F, $_);
	my @a = <F>;
	close F;

	print OUT join( "", @a );
    }
    close OUT;
}


# FIXME: Use qual string and parse base string properly...
sub count_bases{
    my ($b, $q) = @_;
    if ($b) {
	my $a = ($b =~ s/A/A/gi);
	my $t = ($b =~ s/T/T/gi);
	my $g = ($b =~ s/G/G/gi);
	my $c = ($b =~ s/C/C/gi);
	return {'A'=>$a, 'C'=>$c, 'G'=>$g, 'T'=>$t};
    }
    else {
	return {'A'=>0, 'C'=>0, 'G'=>0, 'T'=>0};
    }
}


# Parse and check command line options
sub get_options{
    my %opt = ( 'bed' => $DEFAULT_SNPBED, 'bedxy' => $DEFAULT_XYBED );
    GetOptions( \%opt, 'bam=s', 'bed=s', 'bedxy=s', 'nocheck', 'overwrite', 'out=s', 
		'thread', 'long', 'position' );
    error( "Parameter --bam required", 1, 1 ) unless $opt{bam};
    error( "Bed file $opt{bed} does not exist.", 1 ) if $opt{bed} and !-s $opt{bed};

    unless ($opt{bedxy} eq "none") {
	error( "Sex determination bed file $opt{bedxy} does not exist.", 1 ) if $opt{bedxy} and  !-s $opt{bedxy};
    }
    else {
	$opt{bedxy} = "";
    }

    error( "Parameter --out required.", 1, 1 ) unless $opt{out};

    if (!$opt{overwrite} and ( -s $opt{out}.".genotypes" or -s $opt{out}.".stats" ) ) {
	error( "Output with prefix $opt{out} already exists. Use --overwrite to overwrite.", 1 );
    }
    return %opt;
}


# Display help text
sub display_usage{
    print "provider.pl --bam [METADATA FILE|'FILEMASK'] --out OUT_FILE_PREFIX\n";

    print " REQUIRED\n".
	  "   --bam        Either the path to a metadata file listing bam files or a\n".
	  "                filemask for bam files\n".
	  "   --out        Prefix of output files\n".
	  " OPTIONAL:\n".
	  "   --bed        Path to BED file of SNPs to use for finger printing\n".
	  "                Default: $DEFAULT_SNPBED\n".
	  "   --overwrite  Overwrite any existing files with same file prefix\n".
          "                Default: OFF\n".
	  "   --nocheck    Don't check if files exist or if chromosomes match.\n".
	  "                Default: OFF\n".
	  "   --long       Output genotypes in long format: sampleID[TAB]snpID[TAB]genotype\n".
	  "                Default: OFF\n".
	  "   --position   Output genomic position as output SNP ID instead of BED-defined\n".
	  "                Default: OFF\n".
	  "   --bedxy      Extra bed file for sex determination. Set to 'none' if not desired.\n".
	  "                Default: $DEFAULT_XYBED\n".
	  "   --thread     Run in threaded mode.\n\n";


}

# Get BAM file names from file mask or metadata file.
sub get_bampaths{
    my ($mask, $nocheck) = @_;

    my (@files) = sort glob $mask;
    
    # No file 
    if (@files == 0) {
	error( "No file(s) found: $mask", 1 );	
    }

    # Single file
    elsif (@files == 1) {
	my $fn = $files[0];
	if (-s $fn or $nocheck) {
	    # If bam/sam file, just return file name
	    return ( $fn=>{ 'name'=>$fn } ) if $fn =~ /\.[bs]am$/;

	    # Otherwise assume metadata table. Read and return data.
	    my %files = read_sample_metadata( $fn, $nocheck );
	    return %files;
	}
	else {
	    error( "No file(s) found: $mask", 1 );
	}
    }

    # Multiple files
    else {
	my %files;
	my $cnt = 0;
	foreach (@files) {
	    $cnt ++;
	    $files{$_}->{ name } = $_;
	}
	return %files;
    }
}


# Parse meta data file
sub read_sample_metadata{
    my($fn, $nocheck ) = @_;
    
    my %files;
    my $cnt;

    open( TABLE, $fn );
    while (<TABLE>) {
	chomp;

	next if /^#/ or /^\s*$/;  # Skip comments and empty lines

	my @dat = split /\t/;
	if (-s $dat[0] or $nocheck) {
	    my $name = ($dat[1] or "unknown".++$cnt);
	    $files{ $dat[0] }->{ name } = $name;
	    $files{ $dat[0] }->{ individual } = $name if !$dat[2];
	    $files{ $dat[0] }->{ individual } = $dat[2] if $dat[2];
	    $files{ $dat[0] }->{ sex } = $dat[3] if $dat[3];
	}
	else {
	    error( "ERROR: File not found '$dat[0]'", 1 );
	}
    }
    return %files;
}


# Perform chi2test for Hardy-Weinberg equilibrium
sub HW_chi2test {
    my ($AA, $AB, $BB) = @_;

    return 0 if $AA == 0 or $BB == 0;

    my $N = $AA + $AB + $BB;
    my $p = ($AA*2 + $AB) / (2*$N);
    my $q = ($BB*2 + $AB) / (2*$N);
    my $eAA = ($p ** 2) * $N;
    my $eBB = ($q ** 2) * $N;
    my $eAB = (2 * $p * $q) * $N;
    
    my $chi2 = ($AA-$eAA)**2 / $eAA +
	($BB-$eBB)**2 / $eBB +
	($AB-$eAB)**2 / $eAB;
    
    if ($chi2 < 3.84) {
	return 1;
    }
    return 0;
}


# Check if all the chromosome names in the BED file are present in the BAM
sub matching_chr_names {
    my ($bam, $bed) = @_;
    my (%bam_chr, %bed_chr);

    my @header = `$SAMTOOLS_PATH view -H $bam`;
    foreach (@header) {
	if( /^\@SQ/ ) {
	    my ($name) = ( $_ =~ /\tSN:(.*?)\t/ );
	    $bam_chr{$name} = 1;
	}
    }

    open( BED, $bed );
    while( <BED> ) {
	my @a = split /\t/;
	$bed_chr{$a[0]} = 1;
    }

    foreach my $chr (keys %bed_chr) {
	return 0 unless $bam_chr{ $chr };
    }

    return 1;
}


# Prints error and quits program
sub error{
    my ($msg, $error_code, $print_usage) = @_;
    print STDERR "*** ERROR: $_[0] ***\n\n";
    &display_usage if $print_usage;
    exit $_[1];
}


sub sum{
    my $sum;
    foreach (@_) {
	$sum += $_;
    }
    return $sum;
}


sub large {
    my ($hash, $pos) = @_;

    my $cnt;
    foreach my $i (sort {$hash->{$b}<=>$hash->{$a}} keys %$hash) {
	$cnt++;
	return $i if $cnt == $pos;
    }
    die "Large outside of array";
}
