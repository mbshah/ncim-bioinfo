#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw(uniq);

my $file=$ARGV[0];
my $pairingstatus=$ARGV[1];
my $outdir=$ARGV[2];
my $file_name;
if($file=~/(.*\/)(.*$)/){$file_name=$2;}
#print "$file\t$file_name\n\n";
my $combined_fastq;
#============Combining fastq=========
if ($pairingstatus eq "paired"){
	print "->Joining into a single file\n";
	my $combined_dir=$outdir."/combined_fastq/";
	`mkdir $combined_dir`;
	$combined_fastq=$outdir."/combined_fastq/$file_name"."_combined_fastq.fastq";
	my $r1=$file."_R1.fastq";
	my $r2=$file."_R2.fastq";
	`cat $r1 $r2 >$combined_fastq`;
}
#========No need to combine=======
if ($pairingstatus eq "unpaired"){
	$combined_fastq=$file."\.fastq";
	my $combined_dir=$outdir."/combined_fastq/";
	`mkdir $combined_dir`;
	my $newfile=$outdir."/combined_fastq/$file_name"."_combined_fastq.fastq";
	#print "`cp $combined_fastq $newfile`";
	`cp $combined_fastq $newfile`;
	$combined_fastq=$newfile
}

#========Further Processing is done together========
	print "->Quality and Adapter Trimming\n";
	my $combined_dir=$outdir."/trimgalore_out/";
	`mkdir $combined_dir`;
	`trim_galore -o $outdir/trimgalore_out $combined_fastq >>$outdir/trimgalore_out/report.txt 2>>$outdir/trimgalore_out/report.txt`;
	print "->Converting fastq to fasta\n";
	$combined_dir=$outdir."/combined_fasta/";
	`mkdir $combined_dir`;
	$combined_fastq=~s/\.fastq$/_trimmed\.fq/;
	$combined_fastq=~s/\/combined_fastq\//\/trimgalore_out\//;
	print "$combined_fastq\n";
	my $combined_fasta=$outdir."/combined_fasta/$file_name."."_combined_fasta.fasta";
	`awk '/^@/{gsub(/^@/,">");print;getline;print}' $combined_fastq > $combined_fasta`;
