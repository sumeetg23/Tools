#!/usr/bin/perl
# Created by Sumeet Gupta (sgupta@wi.mit.edu) - 21st June 2016
# Modified by Sumeet Gupta (sgupta@wi.mit.edu) - 6th July 2016
# Modified by Sumeet Gupta (sgupta@wi.mit.edu) - 3rd March 2017 - added demultiplexing task


# Converts FASTQ read header/indentifier format

# Automatically identifies the format of the read id and changes it between the 2 following 2 formats
# read id format 1: @HWUSI-EAS100R:6:73:941:1973#0/1
# read id format 2: @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG

# FASTQ files from WI Genome core end with ";1" or ";2" in format 1 specified above - This script can handle this variation

# Changes
# Version 2
# Added the task option
# Added the ability to generate umi fastq from read 1
# Added the ability to handle paired end reads for format conversion

use warnings;
use strict;
use Carp;
use Storable qw(nstore);
use Cwd;
use Getopt::Long;
use File::Basename;
use Pod::Usage;

my ($task, $read1fastq, $read2fastq,$outputfolder, $help, $options, $man);
my $barcodelength=6;
my $mismatches=2;
my $listofbarcodes = "";
my $randomrunid = int(rand(10)*100);
my $flowcellid = "Unknown";
my $uniquestr = "Processed";

&Getopt::Long::GetOptions(
 		'task=s'		=> \$task,
 		'read1fastq=s'	=> \$read1fastq,
 		'outputfolder=s'	=> \$outputfolder,
 		'uniquestr=s'	=> \$uniquestr,
 		'read2fastq=s'	=> \$read2fastq,
 		'barcodelength=i'	=> \$barcodelength,
 		'flowcellid=s'	=> \$flowcellid,
 		'options'	=> \$options,
  		'man'	=> \$man,
 		'listofbarcodes=s'	=> \$listofbarcodes,
 		'mismatches=s'	=> \$mismatches,
 		'help|h'	=> \$help
) or pod2usage( { -message => "\n\nERROR: Invalid Parameter\n\n" , -verbose => 1} );

if ($help) {
    pod2usage( { -exitstatus => 0 } );
}

if ($options) {
    pod2usage( { -exitstatus => 0 } );
}

if ($man) {
    pod2usage( { -exitstatus => 0, -verbose => 2 } );
}

my $runtype = "SE";

if (!defined $task || (($task ne "changeidformat") && ($task ne "makeread2umi") && ($task ne "demultiplex"))) {
	pod2usage( { -message => "\n\nERROR: Missing or Unrecognized Task\n\n"  , -verbose => 1} );
}
elsif (!defined $read1fastq) {
	pod2usage( { -message => "\n\nERROR: Missing Read 1 File\n\n"  , -verbose => 0} );
}
elsif (!defined $outputfolder) {
	pod2usage( { -message => "\n\nERROR: Output Directory\n\n"  , -verbose => 0} );
}
elsif ($task eq "demultiplex" && $listofbarcodes eq "" ) {
	pod2usage( { -message => "\n\nERROR: List of barcodes separated by a comma is required\n\n"  , -verbose => 0} );
}

if(defined $read2fastq && $read2fastq ne ""){
	$runtype = "PE";
}

my($outfilename, $dirs, $suffix) = fileparse($read1fastq);
$outfilename =~ s/.gz$//g;
$outfilename =~ s/.tar$//g;
$outfilename =~ s/.txt$//g;
$outfilename =~ s/.fq$//g;
my $infilename = $outfilename;
$outfilename = $outputfolder."/".$uniquestr."_".$task."_".$outfilename.".txt";

my $pipecmd = compressioncheck($read1fastq);

if($pipecmd !~ m/cat/){
	print $pipecmd;
	exit;
}

processfq($pipecmd,$task,$randomrunid,$outfilename,$listofbarcodes,$infilename,$outputfolder,$uniquestr,$mismatches);

if($runtype eq "PE") {
	
	($outfilename, $dirs, $suffix) = fileparse($read2fastq);
	$outfilename =~ s/.gz$//g;
	$outfilename =~ s/.tar$//g;
	$outfilename =~ s/.txt$//g;
	$outfilename =~ s/.fq$//g;
	$infilename = $outfilename;
	$outfilename = $outputfolder."/".$uniquestr."_".$task."_".$outfilename.".txt";

	$pipecmd = compressioncheck($read2fastq);
	
	if($pipecmd !~ m/cat/){
		print $pipecmd;
		exit;
	}
	
	processfq($pipecmd,$task,$randomrunid,$outfilename,$listofbarcodes,$infilename,$outputfolder,$uniquestr,$mismatches);
}

sub processfq {
	
	my $pipecmd = shift;
	my $task = shift;
	my $randomrunid = shift;
	my $outputfilename = shift;
	my $listofbarcodes = shift;
	my $infilename = shift;
	my $outputfolder = shift;
	my $uniquestr = shift;
	my $mismatches = shift;
	my @barcodelist=split(/,/,$listofbarcodes);
	
	open(my $PIPEIN, '-|', $pipecmd) or die "Cannot open pipe [$pipecmd]: $!\n";
	
	my %fh;
	my %indexbarcode = ();
	if($task eq "demultiplex") {
			#get all file handles for the barcodes
			my @tempbarcodes = @barcodelist;
			push(@tempbarcodes, 'unknown');
			my $z=0;
			foreach my $outputname  (@tempbarcodes) {
				$tempbarcodes[$z] = $outputname."-".$infilename.".txt";
				$z++;
			}
			%fh = get_write_handles(\@tempbarcodes,$outputfolder,$uniquestr);
	}
	else {
		open(OUTPUT,">$outputfilename") or die "Cannot open output file: $outputfilename: $!\n";
	}
	
	while(!eof($PIPEIN)) {
	
		my ($id, $readseq, $qualityid, $qualityseq, $formatconv, $totalreads) = getnextseq($PIPEIN);
		my ($newid, $barcodeseq) = changeidformat($id,$formatconv,$randomrunid,$flowcellid);
		
		# if no index read is present then the barcodeseq is the first x bases of read 1, x = barcodelength
		if($barcodeseq == 0 || !defined ($barcodeseq)){
			$barcodeseq = substr($readseq,0,$barcodelength);
		}

		my ($newqualityid, $qualityindexeq) = changeidformat($qualityid,$formatconv,$randomrunid,$flowcellid);
		if($task eq "changeidformat"){
			print OUTPUT $newid."\n".$readseq."\n".$newqualityid."\n".$qualityseq."\n";
		}
		elsif($task eq "makeread2umi") {
			my $umiseq = getumiseq($barcodeseq,$barcodelength);
			if($formatconv == 1){
				$id =~ s/\/1\;/\/2\;/g;
				$qualityid =~ s/\/1\;/\/2\;/g;
			}
			elsif($formatconv == 2) {
				$id =~ s/ 1/ 2/g;
				$qualityid =~ s/ 1/ 2/g;
			}
			print OUTPUT $id."\n".$umiseq."\n".$qualityid."\n".substr($qualityseq,0,length($umiseq))."\n";
		}
		elsif($task eq "demultiplex") {
			
			# check if barcode has already been scored
			# if not , then score the barcode
			if(!defined($indexbarcode{$barcodeseq})){
				my %calscore = scorebarcode($barcodeseq,\@barcodelist);
				my $threshold = $barcodelength-$mismatches;
				$indexbarcode{$barcodeseq} = find_barcode_bin(\%calscore,$threshold);
			}
			
			# if yes, then put it in the bin file
			print {$fh{$indexbarcode{$barcodeseq}}} $id."\n".$readseq."\n".$qualityid."\n".$qualityseq."\n";
			
		}
	
	}
	
	close ($PIPEIN);
	close (OUTPUT);
}

sub getumiseq {
	
	my $barcode = shift;
	my $barcodelength = shift;
	
	my $umiseq = "";
	
	if(length($barcode) != $barcodelength){
		$umiseq = substr($barcode,$barcodelength);
	}
	
	return $umiseq;
}

sub getnextseq {
	
	my $PIPEIN =shift;

	my $readid = <$PIPEIN>;
	$readid =~ s/[\n\r]//g;
	
	my $formatconv;
	if($readid =~ m/\#/){
		$formatconv = 1;
	}
	elsif($readid =~ m/\s/) {
		$formatconv = 2;
	}
	
	my $readseq = <$PIPEIN>;
	$readseq =~ s/[\n\r]//g;
	my $qualityid = <$PIPEIN>;
	$qualityid =~ s/[\n\r]//g;
	my $qualityseq = <$PIPEIN>;
	$qualityseq =~ s/[\n\r]//g;
	
	my $totalreads = $.;
	$totalreads = int($totalreads/4);

	return ($readid, $readseq, $qualityid, $qualityseq, $formatconv, $totalreads);
	
}

sub changeidformat {
	
	my ($id,$convtype,$runid,$flowcellid) = @_;
	
	my $formatedid;
	my @seqarray = ();
	my @idstats = ();
	my @barcode = ();
	
	if($convtype == 1){
		
		@seqarray = split(/\#/, $id);
		@barcode = split(/\//, $seqarray[1]);
		@idstats = split(/:/, $seqarray[0]);
		my $readnumber;
		my $passfiler;
		if($barcode[1] =~ m/;/){
			my @readandpf = split(/;/, $barcode[1]);
			$readnumber = $readandpf[0];
			$passfiler ="Y";
			if($readandpf[1] != 1){
				$passfiler = "N";
			}
		}
		else {
			$readnumber = $barcode[1];
			$passfiler ="Y";
		}
		
		$formatedid = $idstats[0].":".$runid.":".$flowcellid.":".$idstats[1].":".$idstats[2].":".$idstats[3].":".$idstats[4]." ".$readnumber.":".$passfiler.":0:".$barcode[0];
		
	}
	elsif($convtype == 2){
		
		@seqarray = split(/\s/, $id);
		@idstats = split(/:/, $seqarray[0]);
		@barcode = split(/:/, $seqarray[1]);
		
		$formatedid = $idstats[0].":".$idstats[3].":".$idstats[4].":".$idstats[5].":".$idstats[6]."#".$barcode[3]."/".$barcode[0];
		
	}
	
	return ($formatedid,$barcode[0]);
	
}

sub compressioncheck {
	
	my ($file) = @_;
	my $pipecmd;
		
	if($file =~ m/.tar.gz$/){
		$pipecmd = "zcat $file | tar  -O -xf -";
	}
	elsif($file =~ m/.gz$/) {
		$pipecmd = "zcat $file"; 
	}
	elsif($file =~ m/.txt$/ || $file =~ m/.fq$/) {
		$pipecmd = "cat $file"; 
	}
	else {
		$pipecmd = "================\n\nFilename does not end with txt or .gz or tar.gz or fq\n\n================\n\nExiting Attempting Format Change!!!!!!!\n\n================\n\n\n";
	}
	
	return $pipecmd;
	
}

sub get_write_handles {
  my @file_names = @{$_[0]};
  my $outputfolder = $_[1];
  my $uniquestr = $_[2];
  my %file_handles;
  foreach (@file_names) {
    open my $fh, '>', $outputfolder."/".$uniquestr."-".$_ or next;
	my @barcordesequence = split('-',$_);
    $file_handles{$barcordesequence[0]} = $fh;
  }
  return %file_handles;
}

sub scorebarcode {
	
	my $barcode = $_[0];
	my @barcodelist = @{$_[1]};
	my %scoring = ();
	foreach my $bar (@barcodelist) {
		my $maxscore = 0;
		my @refbarcode = split(//,$bar);
		my @seqbarcode = split(//,$barcode);
		my $interate = 0;
		foreach my $refbase (@refbarcode) {
			#print "$refbase $seqbarcode[$interate]";<STDIN>;
			if($refbase eq $seqbarcode[$interate]){ 
				$maxscore++;
			}
			$interate++;
		}
		$scoring{$bar}= $maxscore;
	}
	
	return %scoring;
}

sub find_barcode_bin {
	
	my %scoring = %{$_[0]};
	my $threshold = $_[1];
	
	my $barcodeseqtobin;
	my $y=0;
	my @decidescore = ();
	my @blist = ();
	foreach my $key1 (sort { $scoring{$b} <=> $scoring{$a} } keys %scoring) {
		$decidescore[$y] = $scoring{$key1};
		$blist[$y] = $key1;
		$y++;
	}
	
	if($decidescore[0] == $decidescore[1]) {
		$barcodeseqtobin = "unknown";
	}
	else {
		if($decidescore[0] >= $threshold) {  # $threshold = $barcodelength-$barcodemismatches
			$barcodeseqtobin = $blist[0]; 
		}
		else { 
			$barcodeseqtobin = "unknown"; 
		}
	}
	
	return $barcodeseqtobin;
}

__END__

=head1 SYNOPSIS

FASTQ_Format_Conversion_v2.pl [options]

=head1 OPTIONS

=over 4

=item B<--task>

(Required.) User specified task to perform. Input options are "changeidformat" OR "makeread2umi" OR "demultiplex". Reform to script manual (--man option) for details.

=item B<--read1fastq>

(Required.) Full path to the read 1 FASTQ file.

=item B<--outputfolder>

(Required.) Full path to the output folder where the outfile should be written.

=item B<--uniquestr>

(Optional.) Experiment name, which will be used to generate output file names. DEFAULT: "Processed"
    
=item B<--read2fastq>

(Required Only for paired end data) Full path to the read 2 FASTQ file. Required if working with paired end data. Ignore if single end data.

=item B<--barcodelength>

(Optional.) Length of the barcode sequenced. ONLY Required for the "makeread2umi" task. DEFAULT: 6

=item B<--flowcellid>

(Optional.) Flowcell ID on which the sample was sequenced. This is usually present in the run folder name generated by the sequencer OR some FASTQ file read ID formats. Optional. DEFAULT: "Unknown"

=item B<--listofbarcodes>

(Required Only for demultiplex task) List of barcodes separate by a comma. EXAMPLE: TAGTGC,ACGTGA,TGACGT

=item B<--mismatches>

(Optional.) Maximum number of mismatches allowed to map to a barcode. If a sequenced barcode matches to more than 1 expected barcode with the same number of mismatches, the read will not be assigned to either barcode. DEFAULT: 2

=item B<--options>

Prints a brief help message and exits.

=item B<--help>

Print a brief help message and exits.

=item B<--man>

Prints a brief description and help message.

=back

=head1 DESCRIPTION

This program can perform the following 2 tasks:
    
1. "changeidformat" - Converts FASTQ read header/indentifier format

	Automatically identifies the format of the read id and changes it between the 2 following 2 formats
	read id format 1: @HWUSI-EAS100R:6:73:941:1973#0/1
	read id format 2: @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
	
	FASTQ files from WI Genome core end with ";1" or ";2" in format 1 specified above - This script can handle this variation
    
2. "makeread2umi" - Extract the UMI sequence from the index read and put it in a separate FASTQ file. The output can be used by NuGen's tool to identify duplicates using the FASTQ file and alignment file for read 1.

3. "demultiplex" - Demultiplex single FASTQ file into multiple FASTQ files based the expected barcodes to be sequenced.

=cut