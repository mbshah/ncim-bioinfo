#!/usr/bin/env perl
use strict;
use warnings;

#my $ddir="$ARGV[0]/";
my $ddir="/home/ncim/Hot_spring_unk/newprimer/16s FastQ files/merge_out/trimgalore/out/";
opendir my $dir, $ddir or die "Cannot open directory: $ddir $!";
my @in = readdir $dir;
closedir $dir;

open (my $fo,">",$ddir."../randomlyselectedseqs.fasta");
print "picking sequences from the following files: @in";
#@in=("D1_R12.fastq");
foreach my $file(@in){

	if($file=~/fq/){
		print "\now opening $file";
		open (my $ff,"<",$ddir.$file);
		my @fasta;
		my $header;
		my $seq;
		my $header2;
		my $qual;
		while($header=<$ff>){
			$header=~s/@/>/g;
			$seq=<$ff>;
			$header2=<$ff>;
			$qual=<$ff>;
			push (@fasta,$header.$seq);
		}
		print "(".scalar(@fasta)." sequences)\n";
		my @random_set;
		my %seen;

		for (1..20) {
			my $candidate = int rand(scalar(@fasta)-1);
			redo if $seen{$candidate}++;
			push @random_set, $fasta[$candidate];
		}
		my $out= join('', @random_set);
		print $fo $out;
		close($ff);
	}
}
close ($fo);
