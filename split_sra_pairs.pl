#!/usr/bin/env perl
use strict;
use warnings;

my $filepath;
if ( exists($ARGV[0])){
	$filepath=$ARGV[0];
}else{
	
	print "Please Enter File Name:";
	chomp($filepath=<STDIN>);
}


open(my $file, "<",$filepath);
my $fastqlen=`wc -l $filepath`;
$fastqlen=~s/\s.*//;

my $numofseq=$fastqlen/4;
my $flag=1;
$filepath=~s/\.fastq//g;
my $outputR1="$filepath"."_R1.fastq";
my $outputR2="$filepath"."_R2.fastq";

print $filepath;

for (my $i=0;$i<$fastqlen;$i++){
	if($i%4==0 && $i!=0){
		if($flag==1){$flag=2;}
		else{$flag=1;}
		#print "flagis $flag\n";
	}
	my $line=<$file>;
	#print "$flag\n";
	#chomp($line);
	if($flag==1){
		open (my $file, ">>",$outputR1);
		print $file $line;
		close ($file);
	}
	if($flag==2){
		open (my $file, ">>",$outputR2);
		print $file $line;
		close ($file);
	}

}



