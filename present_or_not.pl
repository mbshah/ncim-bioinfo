#!/usr/bin/env perl
#Written by Manan Shah
#contact: m.shah@ncl.res.in or manan.b.shah@gmail.com
use strict;
use warnings;

my %taxtable;

open(my $raw, "<", "/home/ncim/Hot_spring_unk/allsample_swarm/mybiom_no_singletons_use.tsv");
my $header=<$raw>;
$header=<$raw>;
my @sample_names=split("\t",$header);
pop(@sample_names);
shift(@sample_names);
print "@sample_names\n";
my $sample_count=scalar(@sample_names);
print $sample_count;

foreach my $line(<$raw>){
	chomp($line);
	my @line=split("\t",$line);
	shift(@line);
	my $cut_at=scalar(@line)-1;
	my @tax_line=split ("; ",$line[$cut_at]);
	#print "@tax_line\n";
	if (scalar(@tax_line)>=6){
		#print "$tax_line[5]\t";
		my $string="";
		for (my $i=0; $i< scalar(@line)-1;$i++){
			my $symbol;
			if($line[$i]==0){$symbol="-";}else{$symbol="+";}
			$string=$string."$symbol\t";
		}
		#print "$string\n";
		if (exists $taxtable{$tax_line[5]})
		{
			my $table_string=$taxtable{$tax_line[5]};
			my @table_string1=split("\t",$table_string);
			my @new_string=split("\t",$string);
			my $new_string="";
			my $xlim=scalar(@new_string);
			for (my $x=0;$x<$xlim;$x++){
				if($table_string1[$x] ne $new_string[$x]){
					$table_string1[$x]="+";
					
				}
				$new_string=$new_string."$table_string1[$x]\t";
				
			}
			$taxtable{$tax_line[5]}=$new_string;
		}
		else{
		$taxtable{$tax_line[5]}=$string;
		}
	}
	
}
close($raw);
open (my $out, ">", "/home/ncim/Hot_spring_unk/allsample_swarm/present_or_not.tsv");
print $out "Genus\t";
foreach(@sample_names){print $out "$_\t"}
print $out "\n";
foreach(sort(keys %taxtable)){
	print $out "$_\t$taxtable{$_}\n";
}
close($out);
