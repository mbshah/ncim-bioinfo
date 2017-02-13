#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw(uniq);

my $ddir="$ARGV[0]/";
opendir my $dir, $ddir or die "Cannot open directory: $ddir $!";
my @in = readdir $dir;
closedir $dir;
my $outdir="$ARGV[1]merge_out";
`mkdir $outdir`;
`mkdir $outdir/assembled/`;
`mkdir $outdir/unassembled/`;

my @fastqs; my @fastqn; my @fastqj;
foreach(@in){
	if($_=~/fastq/){
		push(@fastqs,$_);
	}
}


foreach(@fastqs){
	if($_=~/(.+)_R.+/){
		push(@fastqn,$1);
	}else{
		my $x=$_;
		$x=~s/.fastq//g;
		push(@fastqj,$x);
	}
}

#print @fastqj;
#print "\n\n";
@fastqn=uniq @fastqn;
open (my $mapfile ,">", "$outdir/../mapping.txt");
print $mapfile "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tReversePrimer\tGroup\tInputFileName\tDescription\n";
foreach(@fastqn){
	
	my $inforward="$ddir"."$_"."_R1.fastq";
	my $inreverse="$ddir"."$_"."_R2.fastq";
	my $log="$outdir/"."$_".".pandaseq.report";
	my $assembled="$outdir"."/assembled/"."$_".".assembled.fasta";
	my $unassembled="$outdir"."/unassembled/"."$_".".unassembled.fasta";
	print "pandaseq -f $inforward  -r $inreverse -g $log -u $unassembled -N -w $assembled\n\n";
	`pandaseq -f $inforward  -r $inreverse -g $log -u $unassembled -N -w $assembled`;
	
	print $mapfile "$_\tTCCCTTGTCTCC\tCCTACGGGNBGCASCAG\tGATCTACNVGGGTATCTAATCC\tA\t$_".".assembled.fasta\t$_\n";
}

print "Converting FastQ to Fasta\n";
foreach(@fastqj){
	my $samplid=$_;
	$samplid=~s/\_/./g;
	my $in="$ddir"."$_".".fastq";
	my $assembled="$outdir"."/assembled/"."$_".".fq2fa.fasta";
	print "`fastq_to_fasta -i $in -o $assembled`\n\n";
	`fastq_to_fasta -i $in -o $assembled`;
	print $mapfile "$samplid\tTCCCTTGTCTCC\tCCTACGGGNBGCASCAG\tGATCTACNVGGGTATCTAATCC\tA\t$_".".fq2fa.fasta\t$_\n";
}
print "Created mapping file mapping.txt\n";

close ($mapfile);
