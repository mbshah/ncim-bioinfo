#!/usr/bin/env perl
use strict;
use warnings;

my $metaphinderdir="/home/ncim/Softwares/MetaPhinder-master/";
my $blastpath="/usr/bin";

#### running metaphinder###
my $inputfasta=$ARGV[0];
my $output_folder=$ARGV[1];

my $metaphinderexec=$metaphinderdir."MetaPhinder.py";
my $metaphinderdb=$metaphinderdir."database/ALL_140821_hr";

`mkdir $output_folder`;

`$metaphinderexec -i $inputfasta -o $output_folder -d $metaphinderdb -b $blastpath`;

###opening metaphinder output###
open(my $meta, "<","$output_folder"."/output.txt")||die "file not found";
my $header=<$meta>;
my %hit_table;
foreach my $line (<$meta>){
	chomp ($line);
	#print "$line\n";
	my @entry=split("\t",$line);
	if ($entry[1] eq "phage"){
		#print @entry;
		my $value= "\|ANI=".$entry[2]."\|Merged_Coverage=".$entry[3]."\|No_of_hits=".$entry[4]."\|seq_len=".$entry[5];
		$hit_table{$entry[0]}=$value;
	}
}
close($meta);
my $keepfile=$output_folder."/virome_contigs.keep";
`echo "\n" >$keepfile`;

foreach (sort keys %hit_table) {
	`echo $_ >>$keepfile`;
}

###Create filtered fasta of virome###
my @key_list=sort keys %hit_table;
open (my $virome_file, ">", $output_folder."/virome_"."$inputfasta");
open (my $original_fasta, "<", $inputfasta)||die "file not found";
while( my $line=<$original_fasta>){
	chomp ($line);
	chomp(my $seq=<$original_fasta>);
	if($line=~/^>(\S+).*/){
		my $contigid=$1;
		if (exists $hit_table{$contigid}){
			print $virome_file ">$contigid $hit_table{$contigid}\n";
			print $virome_file "$seq\n";
		}
	}
}
close ($virome_file);
close ($original_fasta);
