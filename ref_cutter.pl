#!/usr/bin/env perl
use strict;
use warnings;

my $filepath;
my $region;
if ( exists($ARGV[0])){
	$filepath=$ARGV[0];
	print "the refernce file bieng cut:$filepath";
}else{
	
	print "Please enter reference file with full path:\n";
	chomp($filepath=<STDIN>);
	#$filepath="/home/ncim/Softwares/ancillary/Silva_123_release/SILVA123_QIIME_release/rep_set/rep_set_16S_only/97/testref.fasta";
	#$filepath="/home/ncim/Softwares/ancillary/Silva_123_release/SILVA123_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta";
}
if ( exists($ARGV[1])){
	$region=$ARGV[1];
	print "the region bieng cut:$filepath";
}else{
	
	print "Please enter region of interest:\n";
	chomp($region=<STDIN>);
	#$region="V4-V9";
}
$region=uc($region);
my %cut_start=( 
	"V1"=>0,
	"V2"=>100,
	"V3"=>250,
	"V4"=>550,
	"V5"=>750,
	"V6"=>875,
	"V7"=>1075,
	"V8"=>1200,
	"V9"=>1375
);
my %cut_end=( 
	"V1"=>200,
	"V2"=>350,
	"V3"=>625,
	"V4"=>775,
	"V5"=>975,
	"V6"=>1150,
	"V7"=>1275,
	"V8"=>1400,
	"V9"=>1525
);

my @regions_of_interest=split("-",$region);

my $substr_start=$cut_start{$regions_of_interest[0]};
my $substr_end=$cut_end{$regions_of_interest[-1]};
my $substr_len=$substr_end-$substr_start;
if($region eq "XX"){
	print "Cut Start point:";
	chomp($substr_start=<STDIN>);
	print "Cut length:";
	chomp($substr_len=<STDIN>);
	$substr_end=$substr_start+$substr_len;
}
print "Regions needed: @regions_of_interest\nCutting at:$substr_start-$substr_end\nCut Length=$substr_len\n";

my $sum=0;
my $count2=0;
open (my $reffile,"<",$filepath) || die("file not found");
open (my $reffileout,">",$filepath.$region."cut.fasta");
foreach my $line (<$reffile>){
	chomp($line);
	if ($line=~/^>/){
		print $reffileout "$line\n";
		$count2++;
	}
	else{
		my $len=length($line);
		if ($substr_end > $len){$substr_len=$len-$substr_start;}
		my $seq_of_intr=substr($line,$substr_start,$substr_len);
		print $reffileout "$seq_of_intr\n";
		$sum+=$len;
	}
}
my $avg_len=$sum/$count2;
print "Number of Sequences:$count2\nAverage Length of Sequences:$avg_len\n";
close($reffile);
close($reffileout);
