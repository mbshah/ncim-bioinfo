#!/usr/bin/env perl
use strict;
use warnings;
use Statistics::R;

my $infolder;
if ( exists($ARGV[0])){
	$infolder=$ARGV[0];
}else{
	print "enter full path to raw files\n";
	chomp(my $infolder=<STDIN>);#dummy
}
chomp($infolder);
if($infolder=~/\/$/){chop $infolder;}

my $outfolder;
if ( exists($ARGV[1])){
	$outfolder=$ARGV[1];
}else{
	print "enter full path to where you want the output files\n";
	chomp(my $outfolder=<STDIN>);#dummy
}
chomp($outfolder);
if($outfolder=~/\/$/){}else{$outfolder=$outfolder."/";}
my $percent;
if ( exists($ARGV[2])){
	$percent=$ARGV[2];
}else{
	print "Since no Percentage provided, we shall use 97% identity cutoff ";
	chomp(my $percent=97);#dummy
}
my @percent_acceptable=(80,90,94,97,99);
if ( grep( /^$percent$/, @percent_acceptable ) ){}else{print "percent value must be one of \"80,90,94,97,99\"\n";die;}

print "Input Folder:$infolder\nOutput Folder:$outfolder\nPercent Similarity to be matched:$percent\n\n";

print "Enter the full path to Silva Qiime folder:\n";
chomp (my $silva_path="/home/ncim/Softwares/ancillary/Silva_123_release/SILVA123_QIIME_release");#dummy
print "$silva_path\n";

print "Enter the full path to Tax4Fun Silva folder:\n";
chomp (my $tax4fun_silva="/home/ncim/Softwares/ancillary/Tax4Fun_SILVA/SILVA123/");#dummy
print "$tax4fun_silva\n";

`mkdir $outfolder`;
my @args=($infolder,$outfolder);

print "Running\nMerging Paired Ends Read\n";
system ("merge.pl", @args);

my $outfolderx=chop($outfolder);
print "To continue with auto generated mapping file press enter, else replace the mapping file with your own mapping.txt and press enter\n
you may want to change the \"group\" coloumn\n\nsource";
chomp(my $dummy=""); #dummy
print "Preparing for Qiime processing\n";
my $datestring = localtime();
print "\n\n[$datestring]`add_qiime_labels.py -i $outfolder/merge_out/assembled -c InputFileName -o $outfolder/combined -m $outfolder/mapping.txt`";
						`add_qiime_labels.py -i $outfolder/merge_out/assembled -c InputFileName -o $outfolder/combined -m $outfolder/mapping.txt`;



my $repset= $silva_path."/rep_set/rep_set_16S_only/$percent/"."$percent"."_otus_16S.fasta";
my $taxonomy= $silva_path."/taxonomy/16S_only/$percent/raw_taxonomy.txt";
my $similarity=$percent/100;
$datestring = localtime();
print "\n\n[$datestring]Starting QIIME Processing:";
$datestring = localtime();
print "\n\n[$datestring]`pick_otus.py -i $outfolder/combined/combined_seqs.fna -o $outfolder/pick_otus/ -z -s $similarity -m swarm`;";
						`pick_otus.py -i $outfolder/combined/combined_seqs.fna -o $outfolder/pick_otus/ -z -s $similarity -m swarm`;
$datestring = localtime();
print "\n\n[$datestring]`pick_rep_set.py -i $outfolder/pick_otus/combined_seqs_otus.txt -f $outfolder/combined/combined_seqs.fna -r $repset -o $outfolder/pick_otus/st_repset.fna`;";
						`pick_rep_set.py -i $outfolder/pick_otus/combined_seqs_otus.txt -f $outfolder/combined/combined_seqs.fna -r $repset -o $outfolder/pick_otus/st_repset.fna`;
$datestring = localtime();
print "\n\n[$datestring]`assign_taxonomy.py -i $outfolder/pick_otus/st_repset.fna -r $repset -t $taxonomy -o $outfolder/assign_taxonomy/ `;";
						`assign_taxonomy.py -i $outfolder/pick_otus/st_repset.fna -r $repset -t $taxonomy -o $outfolder/assign_taxonomy/ `;
$datestring = localtime();
print "\n\n[$datestring]`make_otu_table.py -i $outfolder/pick_otus/combined_seqs_otus.txt -t $outfolder/assign_taxonomy/st_repset_tax_assignments.txt -o $outfolder/mybiom.biom`;";
						`make_otu_table.py -i $outfolder/pick_otus/combined_seqs_otus.txt -t $outfolder/assign_taxonomy/st_repset_tax_assignments.txt -o $outfolder/mybiom.biom`;
$datestring = localtime();
print "\n\n[$datestring]`biom convert -i $outfolder/mybiom.biom -o $outfolder/mybiom.tsv --to-tsv --header-key taxonomy`;";
						`biom convert -i $outfolder/mybiom.biom -o $outfolder/mybiom.tsv --to-tsv --header-key taxonomy`;
$datestring = localtime();
print "\n\n[$datestring]`parallel_align_seqs_pynast.py -i $outfolder/pick_otus/st_repset.fna -o $outfolder/pynast_aligned_seqs -T --jobs_to_start 2`;";
						`parallel_align_seqs_pynast.py -i $outfolder/pick_otus/st_repset.fna -o $outfolder/pynast_aligned_seqs -T --jobs_to_start 2`;
$datestring = localtime();
print "\n\n[$datestring]`filter_alignment.py -o $outfolder/pynast_aligned_seqs -i $outfolder/pynast_aligned_seqs/st_repset_aligned.fasta`;";
						`filter_alignment.py -o $outfolder/pynast_aligned_seqs -i $outfolder/pynast_aligned_seqs/st_repset_aligned.fasta`;
$datestring = localtime();
print "\n\n[$datestring]`filter_otus_from_otu_table.py -i $outfolder/mybiom.biom -o $outfolder/mybiom_no_singletons.biom -n 2 -s 1`;";
						`filter_otus_from_otu_table.py -i $outfolder/mybiom.biom -o $outfolder/mybiom_no_singletons.biom -n 2 -s 1`;
$datestring = localtime();
print "\n\n[$datestring]`biom convert -i $outfolder/mybiom_no_singletons.biom -o $outfolder/mybiom_no_singletons.tsv --to-tsv --header-key taxonomy`;";
						`biom convert -i $outfolder/mybiom_no_singletons.biom -o $outfolder/mybiom_no_singletons.tsv --to-tsv --header-key taxonomy`;
$datestring = localtime();
print "\n\n[$datestring]`mkdir $outfolder/otus_rarefied`";
						`mkdir $outfolder/otus_rarefied`;
$datestring = localtime();
print "\n\n[$datestring]`single_rarefaction.py -i $outfolder/mybiom_no_singletons.biom -o $outfolder/otus_rarefied/otu_table.biom -d 5000`";
						`single_rarefaction.py -i $outfolder/mybiom_no_singletons.biom -o $outfolder/otus_rarefied/otu_table.biom -d 5000`;
$datestring = localtime();
print "\n\n[$datestring]`filter_fasta.py -f $outfolder/pynast_aligned_seqs/st_repset_aligned_pfiltered.fasta -o $outfolder/pynast_aligned_seqs/st_repset_aligned_biom_filtered_seqs.fastq -b $outfolder/mybiom_no_singletons.biom`;";
						`filter_fasta.py -f $outfolder/pynast_aligned_seqs/st_repset_aligned_pfiltered.fasta -o $outfolder/pynast_aligned_seqs/st_repset_aligned_biom_filtered_seqs.fastq -b $outfolder/mybiom_no_singletons.biom`;
$datestring = localtime();
print "\n\n[$datestring]`make_phylogeny.py -i $outfolder/pynast_aligned_seqs/st_repset_aligned_biom_filtered_seqs.fastq -o $outfolder/rep_set.tre`;";
						`make_phylogeny.py -i $outfolder/pynast_aligned_seqs/st_repset_aligned_biom_filtered_seqs.fastq -o $outfolder/rep_set.tre`;
$datestring = localtime();
print "\n\n[$datestring]`core_diversity_analyses.py -o $outfolder/DiversityAnalyses/ -i $outfolder/otus_rarefied/otu_table.biom -m $outfolder/mapping.txt -t $outfolder/rep_set.tre -e 5000`";
						`core_diversity_analyses.py -o $outfolder/DiversityAnalyses/ -i $outfolder/otus_rarefied/otu_table.biom -m $outfolder/mapping.txt -t $outfolder/rep_set.tre -e 5000`;
$datestring = localtime();
print "\n\n[$datestring] Runing Tax4Fun in R using Script:\n\n";
my $qiimein=$outfolder."/mybiom.tsv";
my$qiimeout=$outfolder."/Tax4FunOutput_mybiom.tsv";
#print $qiimein;

my $R = Statistics::R->new();
my $cmds = <<EOF;
library(Tax4Fun)
myQiimeBiom<-importQIIMEData("$qiimein")
folderReferenceData <- "$tax4fun_silva"
Tax4FunOut <- Tax4Fun(myQiimeBiom, folderReferenceData)
AbundanceProfile <- data.frame(t(Tax4FunOut\$Tax4FunProfile))
write.table(AbundanceProfile,"$qiimeout",sep="\\t",col.names=NA)
print(Tax4FunOut)
EOF

print "$cmds\n\n";
my $out = $R->run($cmds);
$datestring = localtime();
print "\n\n[$datestring]Tax4Fun Output:\n\n$out";
$datestring = localtime();
print "\n\n[$datestring] done\n\n";
