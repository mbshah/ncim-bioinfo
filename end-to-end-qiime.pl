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
print "\n\nadd_qiime_labels.py -i $outfolder/merge_out/assembled -c InputFileName -o $outfolder/combined -m $outfolder/mapping.txt\n\n";
`add_qiime_labels.py -i $outfolder/merge_out/assembled -c InputFileName -o $outfolder/combined -m $outfolder/mapping.txt`;



my $repset= $silva_path."/rep_set/rep_set_16S_only/$percent/"."$percent"."_otus_16S.fasta";
my $taxonomy= $silva_path."/taxonomy/16S_only/$percent/raw_taxonomy.txt";
my $similarity=$percent/100;

print "\n\n`parallel_pick_otus_uclust_ref.py -i $outfolder/combined/combined_seqs.fna -o $outfolder/pick_otus/ -r $repset -z -O 1 -s $similarity`;";
`parallel_pick_otus_uclust_ref.py -i $outfolder/combined/combined_seqs.fna -o $outfolder/pick_otus/ -r $repset -z -O 1 -s $similarity`;

print "\n\n`pick_rep_set.py -i $outfolder/pick_otus/combined_seqs_otus.txt -f $outfolder/combined/combined_seqs.fna -r $repset -o $outfolder/pick_otus/st_repset.fna`;";
`pick_rep_set.py -i $outfolder/pick_otus/combined_seqs_otus.txt -f $outfolder/combined/combined_seqs.fna -r $repset -o $outfolder/pick_otus/st_repset.fna`;

print "\n\n`assign_taxonomy.py -i $outfolder/pick_otus/st_repset.fna -r $repset -t $taxonomy -o $outfolder/assign_taxonomy/ `;";
`assign_taxonomy.py -i $outfolder/pick_otus/st_repset.fna -r $repset -t $taxonomy -o $outfolder/assign_taxonomy/ `;

print "\n\n`make_otu_table.py -i $outfolder/pick_otus/combined_seqs_otus.txt -t $outfolder/assign_taxonomy/st_repset_tax_assignments.txt -o $outfolder/mybiom.biom`;";
`make_otu_table.py -i $outfolder/pick_otus/combined_seqs_otus.txt -t $outfolder/assign_taxonomy/st_repset_tax_assignments.txt -o $outfolder/mybiom.biom`;

print "\n\n`biom convert -i $outfolder/mybiom.biom -o $outfolder/mybiom.tsv --to-tsv --header-key taxonomy`;";
`biom convert -i $outfolder/mybiom.biom -o $outfolder/mybiom.tsv --to-tsv --header-key taxonomy`;

print "\n\n`parallel_align_seqs_pynast.py -i $outfolder/pick_otus/st_repset.fasta -o $outfolder/pynast_aligned_seqs -T --jobs_to_start 2`;";
`parallel_align_seqs_pynast.py -i $outfolder/pick_otus/st_repset.fna -o $outfolder/pynast_aligned_seqs -T --jobs_to_start 1`;

print "\n\n`filter_alignment.py -o $outfolder/pynast_aligned_seqs -i $outfolder/pynast_aligned_seqs/st_repset_aligned.fasta`;";
`filter_alignment.py -o $outfolder/pynast_aligned_seqs -i $outfolder/pynast_aligned_seqs/st_repset_aligned.fasta`;

print "\n\n`make_phylogeny.py -i $outfolder/pynast_aligned_seqs/st_repset_aligned_pfiltered.fasta -o $outfolder/rep_set.tre`;";
`make_phylogeny.py -i $outfolder/pynast_aligned_seqs/st_repset_aligned_pfiltered.fasta -o $outfolder/rep_set.tre`;

print "\n\n`filter_otus_from_otu_table.py -i $outfolder/mybiom.biom -o $outfolder/mybiom_no_singletons.biom -n 2 -s 1`;";
`filter_otus_from_otu_table.py -i $outfolder/mybiom.biom -o $outfolder/mybiom_no_singletons.biom -n 2 -s 1`;

print "\n\n`biom convert -i $outfolder/mybiom_no_singletons.biom -o $outfolder/mybiom_no_singletons.tsv --to-tsv --header-key taxonomy`;";
`biom convert -i $outfolder/mybiom_no_singletons.biom -o $outfolder/mybiom_no_singletons.tsv --to-tsv --header-key taxonomy`;

print "\n\n`mkdir $outfolder/otus_rarefied`";
`mkdir $outfolder/otus_rarefied`;

print "\n\n`single_rarefaction.py -i $outfolder/mybiom_no_singletons.biom -o $outfolder/otus_rarefied/otu_table.biom -d 5000`";
`single_rarefaction.py -i $outfolder/mybiom_no_singletons.biom -o $outfolder/otus_rarefied/otu_table.biom -d 5000`;

print "\n\n`core_diversity_analyses.py -o $outfolder/DiversityAnalyses/ -i $outfolder/otus_rarefied/otu_table.biom -m $outfolder/mapping.txt -t $outfolder/rep_set.tre -e 5000`";
`core_diversity_analyses.py -o $outfolder/DiversityAnalyses/ -i $outfolder/otus_rarefied/otu_table.biom -m $outfolder/mapping.txt -t $outfolder/rep_set.tre -e 5000`;

print "\n\n Runing Tax4Fun in R using Script:\n\n";
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
print "\n\nTax4Fun Output:\n$out";
print "\n\n done\n\n";
