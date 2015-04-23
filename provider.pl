#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;

my %opt = &get_options;

my $SAMTOOLS_PATH = "samtools";
my $SNPBED_PATH = ($opt{bed} or "HPA_1000G_final_38.bed");

my %file_data = &get_bampaths( $opt{bam}, $opt{nocheck} );

my %data = &get_base_freqs_from_bams( \%file_data, $SNPBED_PATH );

&do_genotyping( \%data );

&print_genotype_table( \%data, $opt{out} );



########################

sub print_genotype_table{
    my( $data, $out ) = @_;

    open (GT, ">".$out.".genotypes");

    print GT "loc";
    print GT "\t$_" foreach sort keys %{ ((values %data)[0])->{samples} };
    print GT "\n";
  
    my ($tot_callable, $tot_sites, $HW_pass, $tot_loc);
    foreach my $loc ( sort keys %$data ) {
	print GT $loc;
	foreach my $sid ( sort keys %{ $data->{$loc}->{samples} } ) {
	    my $genotype = $data{$loc}->{samples}->{$sid}->{gt};
	    print GT "\t";
	    print GT defined( $genotype ) ? $genotype : "NA";
	}
	print GT "\n";
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
	
		if ($ap < 0.6) {
		    $gt = join "/", sort( $max,$second );
		}
		elsif ($ap > 0.9) {
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
	    $gt_cnt{ $num{$gtype} }= $gt_cnt{ $gtype };
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
	my $hw_eq = HW_chi2test( $AA, $AB, $BB );
	$data{$loc}->{HW}        = $hw_eq;
	$data{$loc}->{Ncallable} = $AA + $AB + $BB;
	$data{$loc}->{MAF}       = $AB + 2 * $BB;
    }   
}

sub get_base_freqs_from_bams{
    my( $file_data, $snp_fn ) = @_;

    my $pile_out = $opt{out}.".pile.tmp";

    if( $opt{restart} or !-s $pile_out ) {
	print STDERR "Piling up...";
	system( $SAMTOOLS_PATH." mpileup -l ". $snp_fn." ". join( " ", sort keys %$file_data ).
		" > $pile_out 2> $opt{out}.mpileup.log" );
	print STDERR " Done\n";
    }
    else {
	print STDERR "WARNING: Prior pileup file found, resuming with that file. Use --restart to override\n";
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
    my $a = ($b =~ s/A/A/gi);
    my $t = ($b =~ s/T/T/gi);
    my $g = ($b =~ s/G/G/gi);
    my $c = ($b =~ s/C/C/gi);
    return {'A'=>$a, 'C'=>$c, 'G'=>$g, 'T'=>$t};
}


# Parse and check command line options
sub get_options{
    my %opt;
    GetOptions( \%opt, 'bam=s', 'bed=s', 'nocheck', 'restart', 'out=s' );
    error( "Parameter --bam required", 1 ) unless $opt{bam};
    error( "Bed file $opt{bed} does not exist.", 1 ) if $opt{bed} and !-s $opt{bed};
    error( "Parameter --out required.", 1 ) unless $opt{out};

    if (!$opt{restart} and ( -s $opt{out}.".genotypes" or -s $opt{out}.".stats" ) ) {
	error( "Output with prefix $opt{out} already exists. Use --restart to overwrite.", 1 );
    }
    return %opt;
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


# Prints error and quits program
sub error{
    print $_[0]."\n";
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
