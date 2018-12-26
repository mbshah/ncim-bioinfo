#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw(uniq);

print "Running blast file\n";
my $sourcefile; my $outfile;
if (exists $ARGV[0]){
	$sourcefile=$ARGV[0];
}else{
	print "Please enter source file to be used(freqquency):\n";
	chomp($sourcefile="/home/ncim/ganga_sampling/nanopore/run270518/Segregated/amr_out/frequency.txt");
	#chomp($sourcefile=<STDIN>);
}
if (exists $ARGV[1]){
	$outfile=$ARGV[1];
}else{
	print "Please enter output directory to be created:\n";
	#chomp($outfolder=<STDIN>);
	chomp($outfile="/home/ncim/ganga_sampling/nanopore/run270518/Segregated/amr_out/amr_profile.tsv");
}
my @samples;my @family; my @mechanism; my @file;my %family_profile; my %mechanism_profile;
open (my $freq,"<",$sourcefile);
my $column="";
while (my $line=<$freq>){
	chomp($line);
	push (@file, $line);
	my @line=split /\t/,$line;
	
	if($line[0]=~/^#/){
		push(@samples,$line[-1]);
		$column=$line[-1];
		#print "\n$column\n";
	}else{
		my $amr_mech=$line[1]."_".$line[2];
		my $family=$line[1];
		my $read=$line[-1];
		push(@mechanism,$amr_mech);
		push(@family,$line[1]);
		#print "\t$column\t$family\t$read\n";
		if (exists $family_profile{$family}{$column}){
			my $tmp=$family_profile{$family}{$column};
			$tmp=$tmp+$read;
			$family_profile{$family}{$column}=$tmp;
		}else{
			$family_profile{$family}{$column}=$read;
		}
		if (exists $mechanism_profile{$amr_mech}{$column}){
			my $tmp=$mechanism_profile{$amr_mech}{$column};
			$tmp=$tmp+$read;
			$mechanism_profile{$amr_mech}{$column}=$tmp;
		}else{
			$mechanism_profile{$amr_mech}{$column}=$read;
		}
		
	}
}
@family= uniq @family;
@mechanism=uniq @mechanism;
@samples= uniq @samples;

open(my $out,">",$outfile);
print $out "\t";
foreach (@samples){print $out "$_\t";}
print $out "\n";
foreach my $family (@family){
	print $out "$family\t";
	foreach my $sample (@samples){
	my $value=0;
	if (exists $family_profile{$family}{$sample}){
			$value=$value+$family_profile{$family}{$sample};
	}
	print $out "$value\t";
	}
	print $out "\n";
}

print $out "\n\n\n\n\n\t";
foreach (@samples){print $out "$_\t";}
print $out "\n";
foreach my $family (@mechanism){
	print $out "$family\t";
	foreach my $sample (@samples){
	my $value=0;
	if (exists $mechanism_profile{$family}{$sample}){
			$value=$value+$mechanism_profile{$family}{$sample};
	}
	print $out "$value\t";
	}
	print $out "\n";
}

