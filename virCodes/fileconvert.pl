#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw(uniq);

#=========Checking Files in folder=====================
my $ddir="$ARGV[0]/";
#$ddir="/home/ncim/Hot_spring_ind/Meta_Virome/UNK/fastq/";
opendir my $dir, $ddir or die "Cannot open directory: $ddir $!";
my @in = readdir $dir;
closedir $dir;
my $outdir="$ARGV[1]sequence_files";
#$outdir="/home/ncim/Hot_spring_ind/Meta_Virome/UNK/sequence_files";


#===================Sorting files====================
my @fastqs; my @fastqn; my @fastqj;my @fasta;
foreach(@in){
	if($_=~/fastq/){
		push(@fastqs,$_);
	}
	elsif($_=~/fq/){
		push(@fastqs,$_);
	}
	elsif($_=~/fasta/){
		my $fasta=$_;
		$fasta=~s/.fasta//g;
		push(@fasta,$fasta);
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
@fastqn=uniq @fastqn;

#==========Finding location of running code to use in next step===========
use Cwd 'abs_path';
my $code_loc=abs_path($0);
`mkdir $outdir`;
my $run_loc;
if($code_loc=~/(\/.*)(\/.*\.pl$)/){$run_loc=$1;}


#=========Stsrting conversions===================
print "\n----processing FastQ files to Fasta files----\n";
foreach my $file (@fastqj){
	my $file=$ddir.$file;
	my @args=($file,"unpaired",$outdir);
	print "fastamaker.pl- @args\n\n\n";
	system("$run_loc/fastamaker.pl", @args);
}
foreach my $file (@fastqn){
	my $file=$ddir.$file;
	my @args=($file,"paired",$outdir);
	print "fastamaker.pl @args\n\n\n";
	system("$run_loc/fastamaker.pl", @args);

}
