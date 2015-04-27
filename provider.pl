#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;

my $SAMTOOLS_PATH = "samtools";
my $DEFAULT_SNPBED = "HPA_1000G_final_38.bed";
my $MALE_REGION_37 = "Y:2712190";
my $MALE_REGION_38 = "Y:2844077-2844257";

my %opt = &get_options;
my %file_data = &get_bampaths( $opt{bam}, $opt{nocheck} );
my %variant_data = &read_bed( $opt{bed} );

unless ( $opt{ nocheck } or &matching_chr_names( (keys %file_data)[0], $opt{bed} ) ) {
    error( "Chromosome names do not match between BED and BAMs", 1 );
}


my %data = &get_base_freqs_from_bams( \%file_data, $opt{bed} );

&do_genotyping( \%data );

&print_genotype_table( \%data, $opt{out} );
#&find_matching_samples( \%data, \%file_data );


########################

sub find_matching_samples{
    my( $data, $file_data ) = @_;

    my %dist;
    foreach my $samp1 ( keys %$file_data ) {
	my $name1 = $file_data{$samp1}->{name};
	foreach my $samp2 (keys %$file_data) {
	    next if $samp1 eq $samp2;
	    my $name2 = $file_data{$samp2}->{name};
	    unless (defined( $dist{$name2}->{$name1} )) {
		$dist{$name1}->{$name2} = distance( $data, $name1, $name2 );
		print "$name1\t$name2\t".$dist{$name1}->{$name2}."\n";
	    }
	}
    }
}

sub distance{
    my( $data, $id1, $id2 ) = @_;

    my ( $identical, $tot ) = ( 0, 0 );
    foreach my $loc (keys %$data) {
	if( defined($data{$loc}->{samples}->{$id1}->{gt}) and defined($data{$loc}->{samples}->{$id2}->{gt}) ) {
	    $identical++ if $data{$loc}->{samples}->{$id1}->{gt} eq $data{$loc}->{samples}->{$id2}->{gt};
	    $tot++;
	}
    }

    return ($tot-$identical) / $tot;
}

sub read_bed{
    my $bed = shift;
    open( BED, $bed );
    my %var;
    while( <BED> ) {
	chomp;
	my( $chr, $p0, $p1, $id ) = split /\t/;
	$var{$chr.":".$p1} = $id;
    }
    return %var;
}

sub print_genotype_table{
    my( $data, $out ) = @_;

    open (GT, ">".$out.".genotypes");

    if( ! $opt{'long'} ) {
	print GT "loc";
	print GT "\t$_" foreach sort keys %{ ((values %data)[0])->{samples} };
	print GT "\n";
    }
  
    my ($tot_callable, $tot_sites, $HW_pass, $tot_loc);
    foreach my $loc ( sort keys %$data ) {
	my $snp_id = $opt{position} ? $loc : $variant_data{$loc};

	print GT $snp_id if ! $opt{long};

	foreach my $sid ( sort keys %{ $data->{$loc}->{samples} } ) {
	    my $genotype = $data{$loc}->{samples}->{$sid}->{gt};

	    if( $opt{long} ) {
		print GT "$sid\t$sid\t".$snp_id."\t". 
		         plink_gt( $data{$loc}->{samples}->{$sid}->{basecall} ) ."\n";
	    }
	    else {
		print GT "\t". ( defined( $genotype ) ? $genotype : "NA" ) if ! $opt{long};
	    }
	}

	print GT "\n" if ! $opt{long};
	$tot_callable += $data->{$loc}->{Ncallable};
	$tot_sites    += $data->{$loc}->{Ntotal};
	$HW_pass      += $data->{$loc}->{HW};
	$tot_loc      ++;
    }   
    close GT;

    open (STATS, ">".$out.".stats");
    printf STATS "Average callability: %.2f%%\nSites in H-W equilibrium: %d / %d\n", 100*($tot_callable/$tot_sites),
                                                                                     $HW_pass, $tot_loc;
    close STATS;
}


sub plink_gt{
    my $bc = shift;
    return "0\t0" if !defined($bc) or $bc eq "unclear" or $bc eq "lowdata" ;
    return "$bc\t$bc" if length($bc) == 1;
    my ($a1, $a2) = split /\//, $bc;
    return "$a1\t$a1";
}
 

sub do_genotyping{
    my $data = shift;

    foreach my $loc ( keys %$data ) {
	my $num_samples;

	# Do base calling
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
	$data{$loc}->{MAF}       = $AB + 2 * $BB;
    }   
}

sub get_base_freqs_from_bams{
    my( $file_data, $snp_fn ) = @_;

    my $pile_out = $opt{out}.".pile.tmp";

    if( $opt{overwrite} or !-s $pile_out ) {
	print STDERR "Piling up...";
	system( $SAMTOOLS_PATH." mpileup -l ". $snp_fn." ". join( " ", sort keys %$file_data ).
		" > $pile_out 2> $opt{out}.mpileup.log" );
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
	    $frq{$loc}->{samples}->{ $file_data{$id}->{name} }->{bases} = $bases;
	    $frq{$loc}->{samples}->{ $file_data{$id}->{name} }->{depth} = sum( values %{ $bases } );
	}
    }

    return %frq
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
    my %opt = ( 'bed' => $DEFAULT_SNPBED );
    GetOptions( \%opt, 'bam=s', 'bed=s', 'nocheck', 'overwrite', 'out=s', 
		'threads=i', 'long', 'position' );
    error( "Parameter --bam required", 1, 1 ) unless $opt{bam};
    error( "Bed file $opt{bed} does not exist.", 1 ) if $opt{bed} and !-s $opt{bed};
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
	  "   --bam        Either the path to a metadata file listing bam files\n".
	  "                or a filemask for bam files\n".
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
	  "   --position   Output genomic position as identifier instead of rsID\n".
	  "                Default: OFF\n\n";
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
	    $files{ $dat[0] }->{ individual } = $name if $dat[2];
	    $files{ $dat[0] }->{ sex } = $name if $dat[3];
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
