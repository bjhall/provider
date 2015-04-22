#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;

my %opt = &get_options;

my $SAMTOOLS_PATH = "samtools";
my $SNPBED_PATH = ($opt{bed} or "HPA_1000G_final_38.bed");

my %file_data = &get_bampaths( $opt{bam}, $opt{nocheck} );

my %base_freq = &get_base_freqs_from_bams( \%file_data, $SNPBED_PATH );

foreach my $loc ( keys %base_freq ) {
    print $loc;
    foreach my $sample ( sort keys %{$base_freq{$loc}} ) {
	print "\t".
	    ($base_freq{$loc}->{$sample}->{A} or 0).",".
	    ($base_freq{$loc}->{$sample}->{C} or 0).",".
	    ($base_freq{$loc}->{$sample}->{G} or 0).",".
	    ($base_freq{$loc}->{$sample}->{T} or 0);
    }
    print "\n";
}




########################

sub get_base_freqs_from_bams{
    my( $file_data, $snp_fn ) = @_;

    unless (-s "pile.tmp") {
	print STDERR "Piling up...";
	system( $SAMTOOLS_PATH." mpileup -l ". $snp_fn." ". join( " ", sort keys %$file_data )." > pile.tmp 2> /dev/null" );
    }

    print "loc";
    foreach my $fn (sort keys %$file_data ) {
        print "\t".$file_data{$fn}->{name};
    }   
    print "\n";

    
    my %frq;

    open( PILE, "pile.tmp" );
    while (<PILE>) {
	my @p = split /\t/;
	my $loc = $p[0].":".$p[1];
	
	my $i;
	foreach my $id ( sort keys %$file_data ) {
	    $i += 3;
	    my ($depth, $base, $qual) = ( $p[$i], $p[$i+1], $p[$i+2] );
	    my $bases = count_bases( $base, $qual );
	    $frq{$loc}->{ $file_data{$id}->{name} } = $bases;
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
    return {'A'=>$a, 'C'=>$c, 'G'=>$g, 'T'=>$t, 'depth'=>$a+$c+$g+$t};
}

sub get_options{
    my %opt;
    GetOptions( \%opt, 'bam=s', 'bed=s', 'nocheck' );
    error( "Parameter --bam required", 1 ) unless $opt{bam};
    error( "Bed file $opt{bed} does not exist.", 1 ) if $opt{bed} and !-s $opt{bed};
    return %opt;
}


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
	    return ( $fn=>{ 'name'=>'sample1' } ) if $fn =~ /\.[bs]am$/;

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
	    $files{$_}->{ name } = $_;#"sample".$cnt;
	}
	return %files;
    }
}


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
	}
	else {
	    error( "ERROR: File not found '$dat[0]'", 1 );
	}
    }
    return %files;
}


sub error{
    print $_[0]."\n";
    exit $_[1];
}


