#!/usr/bin/env perl
use strict;
use warnings;

print "Running blast file\n";
my $sourcefile; my $outfolder; my $max_algn;my $bcname;my $anotfile;
if (exists $ARGV[0]){
	$sourcefile=$ARGV[0];
}else{
	print "Please enter source file to be used(fasta):\n";
	#chomp($sourcefile="/home/ncim/Downloads/P5556_raw/_sequence_files/15-ks1.bout");
	chomp($sourcefile=<STDIN>);
}
if (exists $ARGV[1]){
	$outfolder=$ARGV[1];
}else{
	print "Please enter output directory to be created:\n";
	#chomp($outfolder=<STDIN>);
	#chomp($outfolder="/home/ncim/Downloads/P5556_raw/_sequence_files/15-ks1_o/files");
}
	
if (exists $ARGV[2]){
	$max_algn=$ARGV[2];
}else{
	print "max possible read length:\n";
	#chomp($outfolder=<STDIN>);
	chomp($max_algn="133.33");
}
my $cutoff=$max_algn-($max_algn*0.25);
#$cutoff=100;
if (exists $ARGV[3]){
	$bcname=$ARGV[3];
}else{
	print "sample name to append:\n";
	chomp($outfolder=<STDIN>);
	#chomp($bcname="15-ks1");
}
if (exists $ARGV[4]){
	$anotfile=$ARGV[4];
}else{
	print "Megares_ancilary folder:\n";
	#chomp($outfolder=<STDIN>);
	chomp($anotfile="/home/ncim/Softwares/ancillary/Megares/megares_annotations_v1.01.csv");
}

print "parameters:
source file:$sourcefile
Output file:$outfolder
Length cutoff:$cutoff";

my $query=$sourcefile;
my $bout=$sourcefile;
$bout=~s/\.fasta$/\.bout/;
$bout=$bout;
$sourcefile=$bout;

`blastn -query $query -db megares.fasta -outfmt '6 qseqid sseqid sacc qcovs evalue pident bitscore qstart qend sstart send length stitle qseq' -out $bout -num_threads 3 -max_target_seqs 1`;
print "\n\n->Loading Anotation file";
my %anotation;
open (my $anot,"<",$anotfile);
while(my $line=<$anot>){
	my @anotvalues=split /,/,$line;
		if($anotvalues[0]=~/Confirmation$/){$anotvalues[2]=$anotvalues[2]."-->requiresSNPconfirmation";}
	$anotation{$anotvalues[0]}="$anotvalues[1]"."_"."$anotvalues[2]";
}

close $anot;

print "
->Creating FASTA Files";

`mkdir $outfolder`;
my %frequncy_table;
my %taxonomy;
open(my $inf, '<',$sourcefile)||warn("file not found");
while(my $line=<$inf>){
	chomp($line);
	my @entry=split /\t/,$line;
	if ($entry[11]>=$cutoff){
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
my %combined_current_annotation;
open (my $freq, ">>", $outfolder."/../frequency.txt");
print $freq "#Sample_name\tAMR_family\tAMR_mechanism\tNumber_of_reads\t$bcname\n";
foreach (sort keys %frequncy_table){
	my $key=$anotation{$taxonomy{$_}};
	my $freq=$frequncy_table{$_};
	if (exists $combined_current_annotation{$key}){
		my $tmp=$combined_current_annotation{$key};
		$tmp=$tmp+$freq;
		$combined_current_annotation{$key}=$tmp;
	}else{
		$combined_current_annotation{$key}=$freq;
	}
	#print $freq "$bcname\t$anotation{$taxonomy{$_}}\t$frequncy_table{$_}\t$taxonomy{$_}\t$_\n";
}
foreach (sort keys %combined_current_annotation){
	my $key=$_;
	$key=~s/_/\t/g;
	print $freq "$bcname\t$key\t$combined_current_annotation{$_}\n";
}
