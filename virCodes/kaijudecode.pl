#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw(uniq);
use experimental 'smartmatch';
use Getopt::Std;

my %options=();
getopts("i:t:p:", \%options);

#DB files
my $taxdump_dir=$options{t} || "/home/manan/NCL/kaiju_decodes/ncbi_tax_dmp/new_taxdump/";
my $fullnamelineagefile=$taxdump_dir."/taxidlineage.dmp";
my $nodesfile=$taxdump_dir."/nodes.dmp";
my $namesfile=$taxdump_dir."/names.dmp";

#input folder
my $input_folder=$options{i};
$input_folder="/home/manan/NCL/kaiju_decodes/upstream/";

#output folder
my $outdir=$options{o} ||"/home/manan/NCL/kaiju_decodes/us4/";

my $lineage_read=$outdir."/read_lineage.tsv";
my $level1_outfile=$outdir."/Domain.tsv";
my $level2_outfile=$outdir."/phylum.tsv";
my $level3_outfile=$outdir."/class.tsv";
my $level4_outfile=$outdir."/order.tsv";
my $level5_outfile=$outdir."/family.tsv";
my $level6_outfile=$outdir."/genus.tsv";
my $classified_outfile=$outdir."/classified.tsv";
my @created_dirs;
my @created_files;

#variable definitions
my %lineage;
my %names;
my %nodes;
my %classified_count;
my %unclassified_count;
my %profile_l1;
my %profile_l2;
my %profile_l3;
my %profile_l4;
my %profile_l5;
my %profile_l6;
my $state;
my %keycheck;
my @barcodes_list;

print "______________________Loading Taxonomy file__________________\n\n";

#lineage
open (my $linfile,"<",$fullnamelineagefile);
foreach my $line(<$linfile>){
	my @entry=split /\|/,$line;
	$entry[0]=~s/\s+$//g;
	$entry[1]=~s/\s+$//g;
	$entry[1]=~s/^\s+//g;
	$lineage{$entry[0]}=$entry[1];
}
close $linfile;
my $db_count=scalar(keys %lineage);
print "Loaded $db_count lineage entries in memory\n\n";
#	foreach (keys %lineage){
#		print "$_\t$lineage{$_}\n";
#	}

#names
$db_count=0;
open (my $namfile,"<",$namesfile);
foreach my $line(<$namfile>){
	my @entry=split /\|/,$line;
	if ($entry[3]=~/scientific name/){
		$entry[0]=~s/\s+$//g;
		$entry[1]=~s/\s+$//g;
		$entry[1]=~s/^\s+//g;
		$names{$entry[0]}=$entry[1];
	}
}
close $namfile;
$db_count=scalar(keys %names);
print "Loaded $db_count names entries in memory\n\n";

#nodes
$db_count=0;
open (my $nodfile,"<",$nodesfile);
foreach my $line(<$nodfile>){
	my @entry=split /\|/,$line;
	$entry[0]=~s/\s+$//g;
	$entry[2]=~s/\s+$//g;
	$entry[2]=~s/^\s+//g;
	$nodes{$entry[0]}=$entry[2];
}
close $nodfile;
$db_count=scalar(keys %nodes);
print "Loaded $db_count nodes entries in memory\n\n";



###################################starting actual processing####################

`mkdir $outdir`;
`mkdir $outdir/segregated`;

my $kaijufiles=`find $input_folder -name *.out`;
my @kaijufiles=split /\n/,$kaijufiles;
open (my $out,">", $lineage_read);
foreach my $kaiju_out(@kaijufiles){
	my $barcode=$kaiju_out;
	$barcode=~s/$input_folder//;
	$barcode=~s/\.out//;
	$keycheck{$barcode}=$barcode;
	if(@barcodes_list~~$barcode){}else{push(@barcodes_list,$barcode);}
	open(my $kout,"<",$kaiju_out);
	while(my $read_result=<$kout>){
		my $lineage;
		my @entries=split /\t/, $read_result;
		chomp(my $classified=$entries[0]);
		chomp(my $read_id=$entries[1]);
		chomp(my $taxid=$entries[2]);
		if ($classified eq "C"){
			if ($taxid == 1){unclassified_counter($barcode);next;}
			if(exists $lineage{$taxid}){$lineage=$lineage{$taxid}." $taxid"}else{$lineage=check_merged($taxid)." $taxid";}
			$lineage=~s/\D/\;/g;
			$lineage=~s/131567//;
			$lineage=~s/^\;//;
			if ($lineage eq ""){unclassified_counter($barcode);next;}else{
				#print "$lineage\t";
				print $out "$read_id\t$barcode\t$lineage\n";
				$lineage=taxonomy_refiner($lineage);
				#print $out "$read_id\t$barcode\t$lineage\n";
				add_to_profiles($barcode,$lineage);
				classified_counter($barcode);
			}
		}else{
			#print "-----";
			unclassified_counter($barcode);
		}
			#print "\n\n\n-n\n\n$read_id\t$taxid\t$lineage";
	}
}



open (my $classi, ">", $classified_outfile);
@barcodes_list= uniq @barcodes_list;

print $classi "barcode\t";
foreach my $bc(sort @barcodes_list){
	print $classi "$bc\t";
}
print $classi "\nsampleID\t";
foreach my $bc(sort @barcodes_list){
	print $classi "$keycheck{$bc}\t";
}
print $classi "\nClassified\t";
foreach my $bc(sort @barcodes_list){
	print $classi "$classified_count{$bc}\t";
}
print $classi "\nUnclassified\t";
foreach my $bc(sort @barcodes_list){
	print $classi "$unclassified_count{$bc}\t";
}
close $classi;



print "\n\n______________________writingclassification_______________\n-->L1";
my $header_all;
my @samples_list=@barcodes_list;
$header_all="sample\t";
open (my $out1,">", $level1_outfile);
print $out1 "sample\t";
foreach my $bc(sort @samples_list){
	print $out1 "$bc\t";
	$header_all=$header_all."$bc\t";
}
print $out1 "\n";
foreach my $org (sort keys %profile_l1){
	print $out1 "$org\t";
	foreach my $bc(sort @samples_list){
		my $profile_l1_score;
		if (exists$profile_l1{$org}{$bc}){$profile_l1_score=$profile_l1{$org}{$bc}}else{$profile_l1_score=0}
		print $out1 "$profile_l1_score\t";
	}
	print $out1 "\n";
}
close $out1;

print "--done\n-->L2";
open (my $out2,">", $level2_outfile);
print $out2 "sample\t";
foreach my $bc(sort @samples_list){
	print $out2 "$bc\t";
}
print $out2 "\n";
foreach my $org (sort keys %profile_l2){
	my $line= "$org\t";
	foreach my $bc(sort @samples_list){
		my $profile_l2_score;
		if (exists$profile_l2{$org}{$bc}){$profile_l2_score=$profile_l2{$org}{$bc}}else{$profile_l2_score=0}
		$line=$line."$profile_l2_score\t";
	}
	level_wise_writer("level_2",$line);
	print $out2 "$line\n";
}
close $out2;

print "--done\n-->L3";
open (my $out3,">", $level3_outfile);
print $out3 "sample\t";
foreach my $bc(sort @samples_list){
	print $out3 "$bc\t";
}
print $out3 "\n";
foreach my $org (sort keys %profile_l3){
	my $line= "$org\t";
	foreach my $bc(sort @samples_list){
		my $profile_l3_score;
		if (exists$profile_l3{$org}{$bc}){$profile_l3_score=$profile_l3{$org}{$bc}}else{$profile_l3_score=0}
		$line=$line."$profile_l3_score\t";
	}
	level_wise_writer("level_3",$line);
	print $out3 "$line\n";
}
close $out3;

print "--done\n-->L4";
open (my $out4,">", $level4_outfile);
print $out4 "sample\t";
foreach my $bc(sort @samples_list){
	print $out4 "$bc\t";
}
print $out4 "\n";
foreach my $org (sort keys %profile_l4){
	my $line= "$org\t";
	foreach my $bc(sort @samples_list){
		my $profile_l4_score;
		if (exists$profile_l4{$org}{$bc}){$profile_l4_score=$profile_l4{$org}{$bc}}else{$profile_l4_score=0}
		$line=$line."$profile_l4_score\t";
	}
	level_wise_writer("level_4",$line);
	print $out4 "$line\n";
}
close $out4;

print "--done\n-->L5";
open (my $out5,">", $level5_outfile);
print $out5 "sample\t";
foreach my $bc(sort @samples_list){
	print $out5 "$bc\t";
}
print $out5 "\n";
foreach my $org (sort keys %profile_l5){
	my $line= "$org\t";
	foreach my $bc(sort @samples_list){
		my $profile_l5_score;
		if (exists$profile_l5{$org}{$bc}){$profile_l5_score=$profile_l5{$org}{$bc}}else{$profile_l5_score=0}
		$line=$line."$profile_l5_score\t";
	}
	level_wise_writer("level_5",$line);
	print $out5 "$line\n";
}
close $out5;

print "--done\n-->L6";
open (my $out6,">", $level6_outfile);
print $out6 "sample\t";
foreach my $bc(sort @samples_list){
	print $out6 "$bc\t";
}
print $out6 "\n";
foreach my $org (sort keys %profile_l6){
	my $line= "$org\t";
	foreach my $bc(sort @samples_list){
		my $profile_l6_score;
		if (exists$profile_l6{$org}{$bc}){$profile_l6_score=$profile_l6{$org}{$bc}}else{$profile_l6_score=0}
		$line=$line."$profile_l6_score\t";
	}
	level_wise_writer("level_6",$line);
	print $out6 "$line\n";
}
close $out6;
print "done\n";





##################################subroutines#####################################
sub taxonomy_refiner{
	my $lineage=$_[0];
	my @lineage=split /\;/,$lineage;
	my $newlineage="";
	if ($names{$lineage[0]} eq "Bacteria"){
		my %lineage;
		my @keys=("superkingdom","phylum","class","order","family","genus","species");
		foreach my $x (@keys){
			$lineage{$x}="-";
		}
		foreach my $i(@lineage){
			if ($nodes{$i}~~@keys){
				$lineage{$nodes{$i}}=$names{$i};
			}
		}
		for (my $y=0;$y<scalar(@keys);$y++){
			if ($lineage{$keys[$y]} eq "-"){$lineage{$keys[$y]}="other_".$lineage{$keys[$y-1]}}
			$newlineage=$newlineage."$lineage{$keys[$y]}; ";
		}
		$newlineage=~s/; $//;
		#print "$newlineage\n\n";
	}
	if ($names{$lineage[0]} eq "Viruses"){
		my %lineage;
		my @keys=("superkingdom","order","family","genus","species");
		foreach my $x (@keys){
			$lineage{$x}="-";
		}
		foreach my $i(@lineage){
			if ($nodes{$i}~~@keys){
				$lineage{$nodes{$i}}=$names{$i};
			}
		}
		for (my $y=0;$y<scalar(@keys);$y++){
			if ($lineage{$keys[$y]} eq "-"){$lineage{$keys[$y]}="other_".$lineage{$keys[$y-1]}}
			$newlineage=$newlineage."$lineage{$keys[$y]}; ";
		}
		$newlineage=~s/; $//;
		#print "$newlineage\n\n";
	}
	if ($names{$lineage[0]} eq "Eukaryota"){
		my %lineage;
		my @keys=("superkingdom","kingdom","phylum","class","order","family","species");
		foreach my $x (@keys){
			$lineage{$x}="-";
		}
		foreach my $i(@lineage){
			if ($nodes{$i}~~@keys){
				$lineage{$nodes{$i}}=$names{$i};
			}
		}
		for (my $y=0;$y<scalar(@keys);$y++){
			if ($lineage{$keys[$y]} eq "-"){$lineage{$keys[$y]}="other_".$lineage{$keys[$y-1]}}
			$newlineage=$newlineage."$lineage{$keys[$y]}; ";
		}
		$newlineage=~s/; $//;
		#print "$newlineage\n\n";
	}
	if ($names{$lineage[0]} eq "Archaea"){
		my %lineage;
		my @keys=("superkingdom","phylum","class","order","family","genus","species");
		foreach my $x (@keys){
			$lineage{$x}="-";
		}
		foreach my $i(@lineage){
			if ($nodes{$i}~~@keys){
				$lineage{$nodes{$i}}=$names{$i};
			}
		}
		for (my $y=0;$y<scalar(@keys);$y++){
			if ($lineage{$keys[$y]} eq "-"){$lineage{$keys[$y]}="other_".$lineage{$keys[$y-1]}}
			$newlineage=$newlineage."$lineage{$keys[$y]}; ";
		}
		$newlineage=~s/; $//;
		#print "$newlineage\n\n";
	}
	
	return $newlineage;
}


sub level_wise_writer{
	my $level=$_[0];
	my $line=$_[1];
	#my $header_all=$_[2];
	#print "$line\n\n";
	if ($line=~/^(\w+);.*/){
		my $out_folder_name=$outdir."/segregated/".$1."/";
		#print "\n$out_folder_name--";
		if ($out_folder_name~~@created_dirs){}else{
			`mkdir $out_folder_name`;
			push (@created_dirs,$out_folder_name);
			@created_dirs=uniq @created_dirs;
		}
		my $outfilename=$out_folder_name.$level.".tsv";
		if($outfilename~~@created_files){}else{
			open (my $file,">",$out_folder_name.$level.".tsv");
			print $file "$header_all\n";
			close $file;
			push(@created_files,$outfilename);
		}
		open (my $file,">>",$out_folder_name.$level.".tsv");
		print $file "$line\n";
		close $file;
	}
}
#print "@created_dirs\n\n\n";

sub classified_counter{
	my $barcode=$_[0];
	if (exists $classified_count{$barcode}){
		my $tmp= $classified_count{$barcode};
		$tmp++;
		$classified_count{$barcode}=$tmp;
	}else{
		$classified_count{$barcode}=1;
	}
}


sub unclassified_counter{
	my $barcode=$_[0];
	if (exists $unclassified_count{$barcode}){
		my $tmp= $unclassified_count{$barcode};
		$tmp++;
		$unclassified_count{$barcode}=$tmp;
	}else{
		$unclassified_count{$barcode}=1;
	}
}

sub check_merged{
	my $taxid=$_[0];
	my $mergedfile=$taxdump_dir."/merged.dmp";
	my $merline= `grep "^$taxid\t" $mergedfile`;
	my $lineage="";
	if ($merline=~/\S+\s+\|\s+(\S+)\s+\|/){$lineage=$lineage{$1}}
	return $lineage;
}


sub add_to_profiles{
	my $barcode=$_[0];
	my $lineage=$_[1];
	my @lineage=split/\;/,$lineage;
	if (exists $lineage[1]){}else{$lineage[1]="other_".$lineage[0]}
	if (exists $lineage[2]){}else{$lineage[2]="other_".$lineage[1]}
	if (exists $lineage[3]){}else{$lineage[3]="other_".$lineage[2]}
	if (exists $lineage[4]){}else{$lineage[4]="other_".$lineage[3]}
	if (exists $lineage[5]){}else{$lineage[5]="other_".$lineage[4]}
	my $level_1=$lineage[0];
	my $level_2=$lineage[0].";".$lineage[1];
	my $level_3=$lineage[0].";".$lineage[1].";".$lineage[2];
	my $level_4=$lineage[0].";".$lineage[1].";".$lineage[2].";".$lineage[3];
	my $level_5=$lineage[0].";".$lineage[1].";".$lineage[2].";".$lineage[3].";".$lineage[4];
	my $level_6=$lineage[0].";".$lineage[1].";".$lineage[2].";".$lineage[3].";".$lineage[4].";".$lineage[5];
	if (exists $profile_l6{$level_6}{$keycheck{$barcode}}){
		my $tmp=$profile_l6{$level_6}{$keycheck{$barcode}};
		$tmp++;
		$profile_l6{$level_6}{$keycheck{$barcode}}=$tmp;
	}else{
		$profile_l6{$level_6}{$keycheck{$barcode}}=1;
	}
	
	if (exists $profile_l5{$level_5}{$keycheck{$barcode}}){
		my $tmp=$profile_l5{$level_5}{$keycheck{$barcode}};
		$tmp++;
		$profile_l5{$level_5}{$keycheck{$barcode}}=$tmp;
	}else{
		$profile_l5{$level_5}{$keycheck{$barcode}}=1;
	}
	
	
	if (exists $profile_l4{$level_4}{$keycheck{$barcode}}){
		my $tmp=$profile_l4{$level_4}{$keycheck{$barcode}};
		$tmp++;
		$profile_l4{$level_4}{$keycheck{$barcode}}=$tmp;
	}else{
		$profile_l4{$level_4}{$keycheck{$barcode}}=1;
	}
	
	
	
	if (exists $profile_l3{$level_3}{$keycheck{$barcode}}){
		my $tmp=$profile_l3{$level_3}{$keycheck{$barcode}};
		$tmp++;
		$profile_l3{$level_3}{$keycheck{$barcode}}=$tmp;
	}else{
		$profile_l3{$level_3}{$keycheck{$barcode}}=1;
	}
	
	
	if (exists $profile_l2{$level_2}{$keycheck{$barcode}}){
		my $tmp=$profile_l2{$level_2}{$keycheck{$barcode}};
		$tmp++;
		$profile_l2{$level_2}{$keycheck{$barcode}}=$tmp;
	}else{
		$profile_l2{$level_2}{$keycheck{$barcode}}=1;
	}
	
	
	if (exists $profile_l1{$level_1}{$keycheck{$barcode}}){
		my $tmp=$profile_l1{$level_1}{$keycheck{$barcode}};
		$tmp++;
		$profile_l1{$level_1}{$keycheck{$barcode}}=$tmp;
	}else{
		$profile_l1{$level_1}{$keycheck{$barcode}}=1;
	}
}
