#!/usr/bin/env perl
#excerpts taken from BLAST eutils manual
BEGIN{$ENV{HTTPS_proxy}="http://172.16.2.37:3128";$ENV{HTTP_proxy}="http://172.16.2.37:3128"}
use strict;
use warnings;
use LWP::Simple;
use List::MoreUtils qw(uniq);


my $db = 'nucleotide';
my $query = '"Viruses"[Organism]+NOT+"cellular+organisms"[Organism]+NOT+wgs[PROP]+NOT+gbdiv+syn[prop]+AND+(srcdb_refseq[PROP]+OR+nuccore+genome+samespecies[Filter])+AND+"complete+genome"[All+Fields]+AND+("3000"[SLEN]+:+"2500000"[SLEN])';

#Search 
my $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
my $url = $base . "esearch.fcgi?db=$db&term=$query&usehistory=y";
print "$url\n";
my $output = get($url);

#parse WebEnv and QueryKey for acquring results
my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
$url = $base . "esummary.fcgi?db=$db&query_key=$key&WebEnv=$web";		#Set base url for further processing


##Acquire Search Summary
open (my $docsum,">", "docsums.txt");
my $retstart=0;
my @taxids;
print "$url\n";
my $retend=0;
while($retend!=1){
my $docsums = get($url."&retstart=$retstart");
print $docsum "$docsums";
print "$retend\t$retstart\n";
$retstart+=10000;
if ($docsums!~/<DocSum>/){$retend=1}
my @temptaxids=($docsums=~m/<Item Name="TaxId" Type="Integer">(\d+)<\/Item>/g);
@temptaxids= uniq @temptaxids;
#print "@temptaxids\n";
if(scalar@temptaxids>=1){@taxids=(@taxids,@temptaxids);}else{$retend=1;}
}
close $docsum;
my $size=scalar(@taxids);
print "taxids: $size";


####Fetching the Fasta Files
#assemble the efetch URL
my $data;
$url = $base . "efetch.fcgi?db=$db&query_key=$key&WebEnv=$web";

#open (my $fasta,">", "av_ml.fasta");
#for(my $ret=0;$ret<=$retstart;$ret+=1000){
#my $url2 = $url."&rettype=fasta&retmode=text&retstart=$ret&retmax=1000";
#print "$ret\t";
#$data = get($url2);
#print $fasta "$data\n";
#}


###Fetching taxonomy lineages
$db='taxonomy';
open (my $tax, ">","taxonomy.tab");
foreach(@taxids){
	my $taxids=$_;
	my $url3="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=$db&id=$taxids&retmode=xml";
	my $taxonomy=get($url3);
	my @lineage=($taxonomy=~/<Lineage>(.*)<\/Lineage>/s);
	my @name=($taxonomy=~/<ScientificName>(.*)<\/ScientificName>/);
	print $tax "$taxids\t$lineage[0]; $name[0]\n";
}
close $tax;
