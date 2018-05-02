#!/usr/bin/env perl
use strict;
use warnings;

my $sourcefile;
if (exists $ARGV[0]){
	$sourcefile=$ARGV[0];
}else{
	print "Please enter source file to be used(blast_ouptufile):\n";
	#chomp($sourcefile="/home/ncim/Student_Data/kushal/ochro_genomes/gbff/OI_LMG_3301.gbk");
	chomp($sourcefile=<STDIN>);
}
my $path;my $file;
if($sourcefile=~/(.*\/)(.+\..+$)/){$path=$1;$file=$2;}else{$path=".";$file=$sourcefile;}
#print "path=$path\nfile=$file\n";

#print "-> creating fasta from gb file\n";

###File operations###
my $gbfilename=$sourcefile;
open (my $gb,"<",$gbfilename);
#my $fastafilename=$gbfilename;
my $fastafilename=$path."../av_fasta.fasta";
#print "reading frm $gbfilename\nwriting to file $fastafilename\n";
open(my $fasta,">>",$fastafilename);


###Variable definition###
my$flag=0;
my $id=0;
my $definition="";
my $sequence="";
my @seq_ids;

    
###Reading File###
while (my $line=<$gb>){
	chomp($line);
	if($flag==1){
		
		if ($line=~/<GBSeqid>(.*)<\/GBSeqid>$/){
			push (@seq_ids,$1);
			if($line=~/<\/GBSeq_other-seqids>/){
				$flag=0;
				$id=join("|",$seq_ids[-1],$seq_ids[0]);
			}
		}
	}
	if($line=~/<GBSeq_other-seqids>/){
		$flag=1;
	}
	if($line=~/<GBSeq_definition>(.*)<\/GBSeq_definition>$/){
		$definition=$1;
	}
	if($line=~/<GBSeq_sequence>(.*)<\/GBSeq_sequence>$/){
		$sequence=$1;
	}
	if($line=~/<\/GBSet>$/){
		$sequence=uc($sequence);
		$id=join("|",$seq_ids[-1],$seq_ids[0]);
		print $fasta ">$id $definition\n$sequence\n";
		my $id=0;
		my $definition="";
		my $sequence="";
	}
}
close $gb;
close $fasta;
#print "\n->done with gb_utils.pl";
