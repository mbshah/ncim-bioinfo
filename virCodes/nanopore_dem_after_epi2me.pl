#!/usr/bin/env perl
use strict;
use warnings;

print "Segregating reads into files\n";
my $sourcefolder; my $outfolder;my %barcode;
if (exists $ARGV[0]){
	$sourcefolder=$ARGV[0];
}else{
	print "Enter EPI2ME downloads FOlder path:\n";
	chomp($sourcefolder="/home/ncim/post_seq/forepi2me/downloads/");
	#chomp($sourcefolder=<STDIN>);
}
if (exists $ARGV[1]){
	$outfolder=$ARGV[1];
}else{
	print "Please enter output directory to be created:\n";
	#chomp($outfolder=<STDIN>);
	chomp($outfolder="/home/ncim/post_seq/run23092018/dem_epi/");
}
my $classified_folder=$sourcefolder."W_OK_CLASS";
my $unclassified_folder=$sourcefolder."W_OK_UNCLASS";
`rm -rf $outfolder`;
`mkdir $outfolder`;

opendir my $dir, $classified_folder or die "Cannot open directory: $classified_folder $!";
my @in_class = readdir $dir;
closedir $dir;

opendir my $dir2, $unclassified_folder or die "Cannot open directory: $unclassified_folder $!";
my @in_unclass = readdir $dir2;
closedir $dir2;

foreach my $file(@in_class){
	if($file eq "."|$file eq ".."){next;}else{
		$file=$classified_folder."/".$file;
		classify($file,$outfolder);	
	}
}
foreach my $file(@in_unclass){
	if($file eq "."|$file eq ".."){next;}else{
		$file=$unclassified_folder."/".$file;
		classify($file,$outfolder);	
	}
}

print "\n\nnumber of sequences for each barcode:\n";
foreach(sort keys %barcode){
	print "$_\t$barcode{$_}\n";
}

print "Creating Fasta:\n";
my $fastafolder=$outfolder."/fastafiles";
`mkdir $fastafolder`;
foreach(sort keys %barcode){
	print "$_";
	my $fastqfolder=$outfolder."/".$_;
	opendir my $dir3, $fastqfolder or die "Cannot open directory: $fastqfolder $!";
	my @in_fastq = readdir $dir3;
	closedir $dir3;
	foreach my $fqf (@in_fastq){
		if($fqf eq "."|$fqf eq ".."){next;}
		my $fastqfile=$fastqfolder."/".$fqf;
		$fqf=~s/\.fastq//;
		my $fastafile=$fastafolder."/".$fqf."_".$_.'.fasta';
		#print "$fastqfile\n$fastafile\n\n";
		`awk '/^@/{gsub(/^@/,">");print;getline;print}' $fastqfile >$fastafile`;
	}
	print "\b\b\b\b\b     \b\b\b\b\b";
}
print "done\n";
=head
print "running blast on fastafiles against All_viruses\n";

#==========Finding location of running code to use in next step===========
use Cwd 'abs_path';
my $code_loc=abs_path($0);
my $run_loc;
if($code_loc=~/(\/.*)(\/.*\.pl$)/){$run_loc=$1;}

#=========================================================================
opendir my $dir4, $fastafolder or die "Cannot open directory: $fastafolder $!";
my @in_fasta = readdir $dir4;
closedir $dir4;
my $boutfolder=$outfolder."/blast_outs/";
`mkdir $boutfolder`;
foreach my $faf (@in_fasta){
	if($faf eq "."|$faf eq ".."){next;}
	my $query=$fastafolder."/".$faf;
	$faf=~s/\.fasta//;
	my $out=$boutfolder.$faf.'.bo';
	my $bpout=$out;
	$bpout=~s/\.bo$/\//;
	my $barcode_bp="";
	if($faf=~/\_(.{2,5})$/){$barcode_bp=$1;}
	print "________________________________________________\nblast $barcode_bp\n";
	`blastn -query $query -db /home/ncim/blastdb/allviruses_ncbi -outfmt '6 qseqid sseqid sacc qcovs evalue pident bitscore qstart qend sstart send length stitle qseq' -out $out -num_threads 3 -max_target_seqs 1`;
	my @args=($out,$bpout,"",$barcode_bp);
	system("$run_loc/blast_parse.pl", @args);
}
=cut

sub classify{
	my $file=$_[0];
	my $outfolder=$_[1];
	print "processing $file to $outfolder\n\n";
	open (my$fastq,"<",$file);
	while(my $line=<$fastq>){
		chomp(my $header=$line);
		chomp(my $sequence=<$fastq>);
		chomp(my $qheader=<$fastq>);
		chomp(my $quality=<$fastq>);
		my $barcode=0;my $runid="";
		if ($header=~/barcode=(\S+)\s*/){$barcode=$1};
		if ($header=~/runid=(\S+)\s*/){$runid=$1};
		#print "$barcode\t$runid\n";
		my $outfile=$outfolder."/".$barcode;
		if(exists $barcode{$barcode}){my $tmp=$barcode{$barcode};$tmp++;$barcode{$barcode}=$tmp;}
		else{`mkdir $outfile`;$barcode{$barcode}=1}		
		$outfile=$outfile."/fastq_runid_".$runid.".fastq";
		open (my $out, ">>", $outfile);
		print $out "$header\n$sequence\n$qheader\n$quality\n";
		close $out;
	}
}
