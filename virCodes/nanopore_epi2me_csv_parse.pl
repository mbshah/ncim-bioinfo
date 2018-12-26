#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw(uniq);
use experimental 'smartmatch';

my $taxdump_dir="/home/ncim/Softwares/ancillary/ncbi_tax_dmp/new_taxdump/";
my $fullnamelineagefile=$taxdump_dir."/fullnamelineage.dmp";
my $csvfile="/home/ncim/ganga_sampling/174714_classification_wimp_v2-v1.csv";
my $outdir="/home/ncim/ganga_sampling/run27072018/epi2me_out2/";
my $epi2me_uploaded="";
my $lineage_read=$outdir."/read_lineage.tsv";
my $extract_list="/home/ncim/ganga_sampling/run27072018/extract.txt";
my $level1_outfile=$outdir."/level1.tsv";
my $level2_outfile=$outdir."/level2.tsv";
my $level3_outfile=$outdir."/level3.tsv";
my $level4_outfile=$outdir."/level4.tsv";
my $level5_outfile=$outdir."/level5.tsv";
my $level6_outfile=$outdir."/level6.tsv";
my $classified_outfile=$outdir."/classified.tsv";
my @created_dirs=undef;

my %lineage;
my %profile_l1;
my %profile_l2;
my %profile_l3;
my %profile_l4;
my %profile_l5;
my %profile_l6;
my @barcodes_list;
my %classified_count;
my %unclassified_count;



print "______________________Loading Taxonomy file in memory__________________\n\n";
open (my $linfile,"<",$fullnamelineagefile);
foreach my $line(<$linfile>){
	my @entry=split /\|/,$line;
	$entry[0]=~s/\s+$//g;
	$entry[2]=~s/\s+$//g;
	$entry[2]=~s/^\s+//g;
	$lineage{$entry[0]}=$entry[2];
}
close $linfile;
my $db_count=scalar(keys %lineage);
print "Loaded $db_count entries in memory\n\n";

my $extractlist=`cat $extract_list`;
my @extract_list=split /\n/, $extractlist;
print "loaded extract list into memory: @extract_list\n";

print "_____________________Processing epi2me csv file___________________________";
`rm -r $outdir`;
`mkdir $outdir`;
open (my $csv,"<", $csvfile);
open (my $out,">", $lineage_read);
my $dummy=<$csv>;
while (my $line=<$csv>){
	my $lineage ="";
	my @entry=split /,/,$line;
	my $barcode=$entry[4];
	if(@barcodes_list~~$barcode){}else{push(@barcodes_list,$barcode);}
	if($entry[3]=~/Classified/){
		my $taxid=$entry[5];
		if(exists $lineage{$taxid}){$lineage=$lineage{$taxid}}else{$lineage=check_merged($taxid);}
		$lineage=~s/cellular organisms;//;
		if ($lineage eq ""){unclassified_counter($barcode);next;}else{
		print $out "$entry[0]\t$barcode\t$lineage\n";
		add_to_profiles($barcode,$lineage);
		classified_counter($barcode);}
	}else{
		unclassified_counter($barcode);
	}
	foreach my $eterm (@extract_list){
		if($entry[6]=~/$eterm/){
			my $orgname=$entry[6];
			$orgname=~s/\s/_/g;
			my $org_dir=$outdir."/".$orgname;
			if (@created_dirs~~$orgname){}else{`mkdir $org_dir`;push(@created_dirs,$orgname);}
			my $outfasta=$org_dir."/".$orgname."_".$barcode.".fastq";
			my $readid=$entry[1];
			my $fastqfilename=$entry[0];
			my $fastqfile=$epi2me_uploaded."/../downloads/W_OK_CLASS/".$fastqfilename;
			#`grep -A3 $readid $fastqfile >>$outfasta`;
		}
	}
}
close $csv;
close $out;

open (my $classi, ">", $classified_outfile);
@barcodes_list= uniq @barcodes_list;
print $classi "\nbarcodes\t";
foreach my $bc(sort @barcodes_list){
	print $classi "$bc\t";
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

print "\n\n______________________writing L1 classification_______________";
open (my $out2,">", $level1_outfile);
print $out2 "\nbarcodes\t";
foreach my $bc(sort @barcodes_list){
	print $out2 "$bc\t";
}
print $out2 "\n";
foreach my $org (sort keys %profile_l1){
	print $out2 "$org\t";
	foreach my $bc(sort @barcodes_list){
		print $out2 "$profile_l1{$org}{$bc}\t";
	}
	print $out2 "\n";
}
close $out2;

print "\n\n______________________writing L2 classification_______________";
open (my $out3,">", $level2_outfile);
print $out3 "\nbarcodes\t";
foreach my $bc(sort @barcodes_list){
	print $out3 "$bc\t";
}
print $out3 "\n";
foreach my $org (sort keys %profile_l2){
	print $out3 "$org\t";
	foreach my $bc(sort @barcodes_list){
		my $profile_l2_score;
		if (exists$profile_l2{$org}{$bc}){$profile_l2_score=$profile_l2{$org}{$bc}}else{$profile_l2_score=0}
		print $out3 "$profile_l2_score\t";
	}
	print $out3 "\n";
}
close $out3;


print "\n\n______________________writing L3 classification_______________";
open (my $out4,">", $level3_outfile);
print $out4 "\nbarcodes\t";
foreach my $bc(sort @barcodes_list){
	print $out4 "$bc\t";
}
print $out4 "\n";
foreach my $org (sort keys %profile_l3){
	print $out4 "$org\t";
	foreach my $bc(sort @barcodes_list){
		my $profile_l3_score;
		if (exists$profile_l3{$org}{$bc}){$profile_l3_score=$profile_l3{$org}{$bc}}else{$profile_l3_score=0}
		print $out4 "$profile_l3_score\t";
	}
	print $out4 "\n";
}
close $out4;

print "\n\n______________________writing L4 classification_______________";
open (my $out5,">", $level4_outfile);
print $out5 "\nbarcodes\t";
foreach my $bc(sort @barcodes_list){
	print $out5 "$bc\t";
}
print $out5 "\n";
foreach my $org (sort keys %profile_l4){
	print $out5 "$org\t";
	foreach my $bc(sort @barcodes_list){
		my $profile_l4_score;
		if (exists$profile_l4{$org}{$bc}){$profile_l4_score=$profile_l4{$org}{$bc}}else{$profile_l4_score=0}
		print $out5 "$profile_l4_score\t";
	}
	print $out5 "\n";
}
close $out5;

print "\n\n______________________writing L5 classification_______________";
open (my $out6,">", $level5_outfile);
print $out6 "\nbarcodes\t";
foreach my $bc(sort @barcodes_list){
	print $out6 "$bc\t";
}
print $out6 "\n";
foreach my $org (sort keys %profile_l5){
	print $out6 "$org\t";
	foreach my $bc(sort @barcodes_list){
		my $profile_l5_score;
		if (exists$profile_l5{$org}{$bc}){$profile_l5_score=$profile_l5{$org}{$bc}}else{$profile_l5_score=0}
		print $out6 "$profile_l5_score\t";
	}
	print $out6 "\n";
}
close $out6;

print "\n\n______________________writing L6 classification_______________";
open (my $out7,">", $level6_outfile);
print $out7 "\nbarcodes\t";
foreach my $bc(sort @barcodes_list){
	print $out7 "$bc\t";
}
print $out7 "\n";
foreach my $org (sort keys %profile_l6){
	print $out7 "$org\t";
	foreach my $bc(sort @barcodes_list){
		my $profile_l6_score;
		if (exists$profile_l6{$org}{$bc}){$profile_l6_score=$profile_l6{$org}{$bc}}else{$profile_l6_score=0}
		print $out7 "$profile_l6_score\t";
	}
	print $out7 "\n";
}
close $out7;





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
			#
	my $barcode=$_[0];
	my $lineage=$_[1];
	my @lineage=split/\;/,$lineage;
	#print $out "$lineage---$lineage[0]";
	if(exists $profile_l1{$lineage[0]}{$barcode}){
		my $tmp=$profile_l1{$lineage[0]}{$barcode};
		$tmp++;
		$profile_l1{$lineage[0]}{$barcode}=$tmp;
	}else{$profile_l1{$lineage[0]}{$barcode}=1}
	
	if (exists $lineage[1]){}else{$lineage[1]="other_".$lineage[0]}
	my $level2=$lineage[0].";".$lineage[1];
	if(exists $profile_l2{$level2}{$barcode}){
		my $tmp=$profile_l2{$level2}{$barcode};
		$tmp++;
		$profile_l2{$level2}{$barcode}=$tmp;
	}else{$profile_l2{$level2}{$barcode}=1}
	
	if (exists $lineage[2]){}else{$lineage[2]="other_".$lineage[1]}
	my $level3=$lineage[0].";".$lineage[1].";".$lineage[2];
	if(exists $profile_l3{$level3}{$barcode}){
		my $tmp=$profile_l3{$level3}{$barcode};
		$tmp++;
		$profile_l3{$level3}{$barcode}=$tmp;
	}else{$profile_l3{$level3}{$barcode}=1}
	
	if (exists $lineage[3]){}else{$lineage[3]="other_".$lineage[2]}
	my $level4=$lineage[0].";".$lineage[1].";".$lineage[2].";".$lineage[3];
	if(exists $profile_l4{$level4}{$barcode}){
		my $tmp=$profile_l4{$level4}{$barcode};
		$tmp++;
		$profile_l4{$level4}{$barcode}=$tmp;
	}else{$profile_l4{$level4}{$barcode}=1}
	
	if (exists $lineage[4]){}else{$lineage[4]="other_".$lineage[3]}
	my $level5=$lineage[0].";".$lineage[1].";".$lineage[2].";".$lineage[3].";".$lineage[4];
	if(exists $profile_l5{$level5}{$barcode}){
		my $tmp=$profile_l5{$level5}{$barcode};
		$tmp++;
		$profile_l5{$level5}{$barcode}=$tmp;
	}else{$profile_l5{$level5}{$barcode}=1}
	
	if (exists $lineage[5]){}else{$lineage[5]="other_".$lineage[4]}
	my $level6=$lineage[0].";".$lineage[1].";".$lineage[2].";".$lineage[3].";".$lineage[4].";".$lineage[5];
	if(exists $profile_l6{$level6}{$barcode}){
		my $tmp=$profile_l6{$level6}{$barcode};
		$tmp++;
		$profile_l6{$level6}{$barcode}=$tmp;
	}else{$profile_l6{$level6}{$barcode}=1}
}

