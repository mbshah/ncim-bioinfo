#!/usr/bin/env perl
use strict;
use warnings;

print "Making Fasta ready for making DB\n";
my $sourcefile; my $outfile;
if (exists $ARGV[0]){
	$sourcefile=$ARGV[0];
}else{
	print "Please enter source file to be used(fasta):\n";
	chomp($sourcefile=<STDIN>);
}
if (exists $ARGV[1]){
	$outfile=$ARGV[1];
}else{
	print "Please enter output file to be created(fasta):\n";
	chomp($outfile=<STDIN>);
}
my $outfolder;
if($outfile=~/(\/.*\/)(.*$)/){$outfolder=$1;}else{$outfolder="./";}
print "$outfile\n$outfolder";

print "parameters:
source file:$sourcefile
Output file:$outfile
Output folder: $outfolder


->Starting filteration\n";

open(my $inf, '<',$sourcefile);
open(my $out, '>',$outfile);
while(my $header=<$inf>){
	chomp($header);
	if($header!~/^>/){
		die "requires single line fasta file which can be made using the command:
		awk \'/^>/ {printf(\"\\n%s\\n\",\$0);next; } { printf(\"%s\",\$0);}\' < all_virus_sequences.fasta |grep \"^\\S\" >all_virus_sequences_singleline.fasta\nexiting";
	}
	my $seq=<$inf>;
	chomp($seq);
	if ($header=~/complete genome/){
		print $out "$header";
		#print "$header\n";
		if($header=~/phage/){
			print $out "------phage";
		}
		my $len=length($seq);
		print $out "------$len\n$seq\n";
	}
}
close $inf;
close $out;
print "->Finished filteration\n";

print "->Making blastdb\n";
my $dbname=$outfolder."allviruses_ncbi";
`makeblastdb -in $outfile -input_type fasta -dbtype nucl -out $dbname -parse_seqids`;

print "->done with dbmaker.pl";
