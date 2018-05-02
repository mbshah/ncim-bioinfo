#!/usr/bin/env perl
use strict;
use warnings;

print "Parsing blast output file\n";
my $sourcefile; my $outfolder; my $max_algn;
if (exists $ARGV[0]){
	$sourcefile=$ARGV[0];
}else{
	print "Please enter source file to be used(blast_ouptufile):\n";
	chomp($sourcefile="/home/ncim/Hot_spring_ind/Meta_Virome/UNK/myblast_out2");
	#chomp($sourcefile=<STDIN>);
}
if (exists $ARGV[1]){
	$outfolder=$ARGV[1];
}else{
	print "Please enter output directory to be created:\n";
	#chomp($outfolder=<STDIN>);
	chomp($outfolder="/home/ncim/Hot_spring_ind/Meta_Virome/UNK/viral_fastqs");
}
if (exists $ARGV[2]){
	$max_algn=$ARGV[2];
}else{
	print "max possible read length:\n";
	#chomp($outfolder=<STDIN>);
	chomp($max_algn="100");
}
my $cutoff=$max_algn-($max_algn*0.25);
print "parameters:
source file:$sourcefile
Output file:$outfolder
Length cutoff:$cutoff

->Creating FASTA Files";

`rm -rf $outfolder`;
`mkdir $outfolder`;
my %frequncy_table;
my %taxonomy;
open(my $inf, '<',$sourcefile)||warn("file not found");
while(my $line=<$inf>){
	chomp($line);
	my @entry=split /\t/,$line;
		if ($entry[11]>=75){
		my $entry_filename=$outfolder."/".$entry[2]."\.fasta";
		my $definition=">".$entry[0]."---".$entry[3]."---".$entry[4]."---".$entry[5]."---".$entry[6]."---".$entry[7]."---".$entry[8];
		my $seq=$entry[13];
		open(my $fasta,">>",$entry_filename);
		print $fasta "$definition\n$seq\n";
		close $fasta;
		if (exists $frequncy_table{$entry[2]}){
			my $temp=$frequncy_table{$entry[2]};
			$temp++;
			$frequncy_table{$entry[2]}=$temp;
		}else{
			$frequncy_table{$entry[2]}=1;
			$taxonomy{$entry[2]}=$entry[12];
		}
	}
}
close $inf;
print "Created \n";
open (my $freq, ">", $outfolder."/../frequency.txt");
foreach (sort keys %frequncy_table){
	print $freq"$taxonomy{$_}\t$frequncy_table{$_}\n";
}
