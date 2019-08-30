#!/usr/bin/env perl
#Written by Manan Shah v2
#contact: m.shah@ncl.res.in or manan.b.shah@gmail.com
#Add comments here to list changelog
#excerpts taken from BLAST eutils manual

#BEGIN{$ENV{HTTPS_proxy}="http://172.16.2.37:3128";$ENV{HTTP_proxy}="http://172.16.2.37:3128"}
use strict;
use warnings;
use LWP::Simple;
use List::MoreUtils qw(uniq);

my $filter='yes';
my $ddir=".";
if (exists $ARGV[0]){
	$ddir=$ARGV[0];
}else{
	print "Please enter outputfolder to be used:\n\n";
	chomp($ddir="/home/manan/NCL/test");
	#chomp($ddir=<STDIN>);
}
if ($ddir!~/\/$/){$ddir=$ddir."/";}
mkdir $ddir;

use Cwd 'abs_path';
my $code_loc=abs_path($0);
my $run_loc;
if($code_loc=~/(\/.*)(\/.*\.pl$)/){$run_loc=$1;}

my $db = 'nucleotide';
my $query = '"Viruses"[Organism]+NOT+"cellular+organisms"[Organism]+NOT+wgs[PROP]+NOT+gbdiv+syn[prop]+AND+(srcdb_refseq[PROP]+OR+nuccore+genome+samespecies[Filter])+AND+"complete+genome"[All+Fields]+AND+("3000"[SLEN]+:+"2500000"[SLEN])';
#$query='"viruses"[porgn]+NOT+partial[All+Fields]+AND+("3000"[SLEN]+:+"2500000"[SLEN])';
#Search 
print "->Searching NCBI Using URL:";
my $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
my $url = $base . "esearch.fcgi?db=$db&term=$query&usehistory=y";
print "$url\n\n";
my $output = get($url);
#print $output;
open (my $searchres,">",$ddir."searchresult.txt");
print $searchres $output;
close $searchres;

#parse WebEnv and QueryKey for acquring results
my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/) or die;
my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
my $rescount=$1 if ($output =~ /<eSearchResult><Count>(\d+)<\/Count>/);
$url = $base . "esummary.fcgi?db=$db&query_key=$key&WebEnv=$web";		#Set base url for further processing


##Acquire Search Summary
print "\n\n->retreiving DocSums from NCBI\n\n";
my $docsumfile=$ddir."docsums.txt";
chomp(my $docsumscount=`grep "<DocSum>" $docsumfile |wc -l`);
#print "$docsumscount\t$rescount\n\n";
if ($docsumscount!=$rescount){

	open (my $docsum,">", $docsumfile);
	my $retstart=0;
	my @taxids;
	print "$url\n\n";
	my $retend=0;
	while($retend!=1){
	my $docsums = get($url."&retstart=$retstart");
	print $docsum "$docsums";
	print "$retend\t$retstart\n\n";
	$retstart+=10000;
	if ($docsums!~/<DocSum>/){$retend=1}
	my @temptaxids=($docsums=~m/<Item Name="TaxId" Type="Integer">(\d+)<\/Item>/g);
	@temptaxids= uniq @temptaxids;
	#print "@temptaxids\n\n";
	if(scalar@temptaxids>=1){@taxids=(@taxids,@temptaxids);}else{$retend=1;}
	}
	close $docsum;
}else {print "---->docsums exist\n\n"}

####Processing Docsumsfile###
print "->Processing DocSum file\n\n";
open (my $ds_r,"<", $ddir."docsums.txt");
my %virus_table;
my $org_id;
my $taxid;
my $length;

print "--->filtering non unique taxa IDs\n\n";
open(my $duplicates ,">", $ddir."duplicates.txt");						#to record duplicates
while(my $line=<$ds_r>){
	if ($line=~/<Item Name="Caption" Type="String">(\w\w?\w?\d+)<\/Item>/){$org_id=$1;}
	if ($line=~/<Item Name="TaxId" Type="Integer">(\d+)<\/Item>/){$taxid=$1;}
	if ($line=~/<Item Name="Length" Type="Integer">(\d+)<\/Item>/){$length=$1;}
	if ($line=~/^<\/DocSum>$/){
		if (exists $virus_table{$taxid}){
			my $temp=$virus_table{$taxid};
			$temp=$temp."; $org_id";
			$virus_table{$taxid}=$temp;
			print $duplicates "$taxid\t$org_id\t$length\n\n";
		}else{
			$virus_table{$taxid}=$org_id;
		}
		$org_id="";
		$taxid="";
	}
}
close $duplicates;


####Preparing to download and format the files
print "--->preparing unique taxa IDs for retrival from NCBI\n\n";
open (my $uniques ,">",$ddir."tax_acc.tab");
my @toret;
foreach (sort keys %virus_table){
	my $temp=$virus_table{$_};
	print $uniques "$_\t$temp\n\n";
	my @orgs=split (';\s',$temp);
	if ($filter eq 'yes'){
		push(@toret,$orgs[0]);
	}else{
		push (@toret,@orgs);
	}
}
close $uniques;


##Making and checking of the files already exist
my $db_folder=$ddir."db/";
`mkdir $db_folder`;
opendir my $dir, $db_folder or die "Cannot open directory: $ddir $!";
my @check = readdir $dir;
closedir $dir;
my @toret2;
foreach(@toret){
	my $check=$_.".gb.xml";
	if ($check~~@check){
		print "skipping $check\n\n";
	}else{
		push(@toret2,$_);
	}
}

my $al=scalar(@toret2);
print "->Now Retriving GB files from NCBI\n\n";
print "--->Number of entries to fetch:$al\n\n";
my $status=retrive_data(\@toret2,$al,$db_folder);						##Actual retreiving step
if ($status eq "done"){print "--->done acquiring\n\n";}else{print "--->failed\n\n";}


###Creating fasta files for making blast db###
print "->Creating FASTA files\n\n";
my $fastafile=$ddir."av_fasta.fasta";
open(my $fasta,">",$fastafile);
close $fasta;
foreach(@toret){
	my $file=$db_folder.$_.".gb.xml";
	my @utilsargs=($file);
	system("$run_loc/gb_utils_xml.pl",@utilsargs);
}

my $dboutfile=$ddir."blastdb/av_fasta_filtered.fasta";
my @dbargs=($fastafile,$dboutfile);
system ("$run_loc/dbmaker.pl",@dbargs);

#subroutine to retreive gb xml files from NCBI
sub retrive_data {
	my $orgs=$_[0];
	my $length=$_[1];
	my $dbf=$_[2];
	my @rerun;
	my $status;
	#print $orgs->[2];
	print "$length\tfetching\n\n";
	for (my $i=0;$i<$length;$i++){
		my $acc=$orgs->[$i];
		#print $acc;
		my $url3="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=$acc&retmode=xml";
		my $entry=get($url3);
		if ($entry=~/<\/GBSet>/){
			my $filename=$dbf.$acc.".gb.xml";
			open (my $record,">",$filename);
			print $record $entry;
			close $record;
		}else{
			push(@rerun,$acc);
		}
		
	}
	if (scalar(@rerun)>0){
			$status=retrive_data(\@toret,$al,$dbf);
			if ($status eq "done"){return "done";}else{return "failed";}
		}else{
			return "done";
		}
	
}
