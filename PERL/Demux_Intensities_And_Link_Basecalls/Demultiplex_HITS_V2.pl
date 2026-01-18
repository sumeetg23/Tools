#!/usr/bin/perl

#######################################################################
# Created by Sumeet Gupta (sgupta@wi.mit.edu) - Dec 2010
# Modified to handle zip files - April 2020
# Input file - assumed 4th column is the barcode in a tab delimited file
#######################################################################

use warnings;
use strict;
use Carp;
use File::Basename;


my $USAGE = "\nperl Demultiplex_HITS_V1.pl <Input file> <LIST OF BARCODES SEPARATED BY COMMA> <NUMBER OF MISMATCHES ALLOWED IN BARCODE> <Output Directory>\n\nInput file - assumed 4th column is the barcode in a tab delimited file\n";
unless(@ARGV == 4) { croak "$USAGE\n"; exit; }

my($ofile, $dirs, $suffix) = fileparse($ARGV[0]); # Incorporated for maybe change in future to have a better way to get the arguments and have output directory as optional
$ARGV[3] =~ s/[\n\r]//g;
$dirs = $ARGV[3]; 

# Open input file
$ARGV[0] =~ s/[\n\r]//g;

if ($ARGV[0] =~ /.gz$/) {
	open(SEQUENCEFILE, "gunzip -c $ARGV[0] |") || die "can’t open pipe to $ARGV[0]";
	$ofile =~ s/.gz$//g;
}
else {
	open(SEQUENCEFILE, "$ARGV[0]") || die "can’t open $ARGV[0]";
}

# paramters
my @barcodelist = split(/,/,$ARGV[1]);
my %indexbarcode = ();
my $barcodelength = length($barcodelist[0]);
my $barcodemismatches = $ARGV[2];

# generate output file handles

my @tempbarcodes = @barcodelist;
push(@tempbarcodes, 'unknown');
my $z=0;
foreach my $outputname  (@tempbarcodes) {
	$tempbarcodes[$z] = $dirs."/".$outputname."-".$ofile;
	print $tempbarcodes[$z];
	$z++;
}
my %fh = get_write_handles(@tempbarcodes);

# Demultiplex
my $barcodepos = 3; #offset by 1 because the index starts from 0

while(my $newline = <SEQUENCEFILE>) {

	$newline =~ s/[\n\r]//g;
	if($newline eq "") {next;} # skip any empty lines - usually only present in custom format files
	my @seqarray = split(/\t/, $newline);
	$seqarray[$barcodepos] =~ s/\./N/g; # converting any . into N's in the barcode.
	if (scalar(@seqarray) > 1) {
	
		my $barcode = substr($seqarray[$barcodepos],0,$barcodelength);

		if(!defined($indexbarcode{$barcode})){
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
			my $y=0;
			my @decidescore = ();
			my @blist = ();
			foreach my $key1 (sort { $scoring{$b} <=> $scoring{$a} } keys %scoring) {
				$decidescore[$y] = $scoring{$key1};
				$blist[$y] = $key1;
				
				$y++;
			}
			
			if($decidescore[0] == $decidescore[1]) {
				$indexbarcode{$barcode} = "unknown";
			}
			else {
				if($decidescore[0] >= $barcodelength-$barcodemismatches) { $indexbarcode{$barcode} = $blist[0]; }
				else { $indexbarcode{$barcode} = "unknown"; }
			}
		}

		print {$fh{$indexbarcode{$barcode}}} $newline."\n";

	}

}

close (SEQUENCEFILE);

# close all file handles
if (scalar(keys(%fh)) > 1 ) {
	foreach my $key (keys %fh) {
		close($fh{$key});
	}
}
else {
	print "ERROR - NO COUNTS";
}

#just output a .finished file to indicate the process completed
my $outputfinish = $ARGV[0]."finished";
open(FINISH,">$outputfinish") or die "Cannot open output file: $outputfinish: $!\n";
close(FINISH);

sub get_write_handles {
  my @file_names = @_;
  my %file_handles;
  foreach (@file_names) {
    open my $fh, '>', $_ or next;
	my @barcordesequence = split('-',basename($_));
    $file_handles{$barcordesequence[0]} = $fh;
  }
  return %file_handles;
}
