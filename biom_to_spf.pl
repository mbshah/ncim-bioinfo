#!/usr/bin/env perl
use strict;
use warnings;


my $infile;
my $outfile;


if ( exists($ARGV[0])){
	$infile=$ARGV[0];
}else{
	print "enter full path to biom tsv files\n";
	#chomp($infile=<STDIN>);#dummy
	chomp($infile="/home/ncim/Hot_spring_unk/T6_OUT2/Taxa_plots/mybiom_no_singletons_use_L50.txt");#dummy
}

if ( exists($ARGV[1])){
	$outfile=$ARGV[1];
}else{
	print "enter full path to outfile\n";
	#chomp($outfile=<STDIN>);#dummy
	chomp($outfile="/home/ncim/Hot_spring_unk/T6_OUT2/Taxa_plots/mybiom_no_singletons_use_L50.spf");#dummy
}

open ( my $biom,"<",$infile)||die "file not found";
open ( my $spf,">",$outfile);
my $flag=3;
my$header="";
my @entry;
my @levels;
my @samples;
my $uncultured_count=1;
my $unknown=1;
my %hierarchy_check;
my %hierarchy_change;
my $autonumber=1;
while (my $line=<$biom>){
	my $counts_line="";
	chomp $line;
	 @entry=split /\t/,$line;
	
	#print "$entry[0]\n";
	@levels=split/;/,$entry[0];
	if ($flag!=0){
		if ($flag==2){@samples=@entry;}
		if($flag==1){
		my $levels=scalar(@levels);
		print "Levels = $levels\n\n\n";
		for(my $i=0;$i<$levels;$i++){
			$header=$header."Level_$i\t";
		}
		shift(@samples);
		foreach(@samples){
			$header=$header.$_."\t";
		}
		chop ($header);
		print $spf $header."\n";
	}
		$flag--;
	}
	for (my $j=0;$j<scalar(@levels);$j++){
		my $val=$levels[$j];
		if($val eq "Other"){
			$levels[$j]="$val $levels[$j-1]";
			my $val=$levels[$j];
			if($j>=2 && exists $hierarchy_check{$val} && $hierarchy_check{$val} ne $levels[$j-1]){
				print "error $val<=>$levels[$j-1] doesnot match $val<=>$hierarchy_check{$val} replacing with";
				my $tochange="$val;$levels[$j-1]";
				if(exists $hierarchy_change{$tochange}){
					$val=$hierarchy_change{$tochange};
					
				}else{
					$hierarchy_change{$tochange}= "$val $autonumber";
					$val=$hierarchy_change{$tochange};
					$autonumber++;
				}
				print " $val--(other)\n";
			}else{
				$hierarchy_check{$val} = $levels[$j-1];
			}
			$counts_line=$counts_line."$val\t";
		}elsif($val eq "uncultured bacterium"){
			$levels[$j]="$val $uncultured_count";
			my $val=$levels[$j];
			
			$uncultured_count++;
			if($j>=2 && exists $hierarchy_check{$val} && $hierarchy_check{$val} ne $levels[$j-1]){
				print "error $val<=>$levels[$j-1] doesnot match $val<=>$hierarchy_check{$val} replacing with";
				my $tochange="$val;$levels[$j-1]";
				if(exists $hierarchy_change{$tochange}){
					$val=$hierarchy_change{$tochange};
					
				}else{
					$hierarchy_change{$tochange}= "$val $autonumber";
					$val=$hierarchy_change{$tochange};
					$autonumber++;
				}
				print " $val--(U bact)\n";
			}else{
				$hierarchy_check{$val} = $levels[$j-1];
			}
			$counts_line=$counts_line."$val\t";
		}elsif($val eq "uncultured"){
			$levels[$j]="$val $uncultured_count";
			my $val=$levels[$j];
			
			$uncultured_count++;
			if($j>=2 && exists $hierarchy_check{$val} && $hierarchy_check{$val} ne $levels[$j-1]){
				print "error $val<=>$levels[$j-1] doesnot match $val<=>$hierarchy_check{$val} replacing with";
				my $tochange="$val;$levels[$j-1]";
				if(exists $hierarchy_change{$tochange}){
					$val=$hierarchy_change{$tochange};
					
				}else{
					$hierarchy_change{$tochange}= "$val $autonumber";
					$val=$hierarchy_change{$tochange};
					$autonumber++;
				}
				print " $val --(uncultured)\n";
			}else{
				$hierarchy_check{$val} = $levels[$j-1];
			}
			$counts_line=$counts_line."$val\t";
		}elsif($val=~/Unknown.*/){
			$levels[$j]="$val $unknown";
			my $val=$levels[$j];
			
			$unknown++;
			if($j>=2 && exists $hierarchy_check{$val} && $hierarchy_check{$val} ne $levels[$j-1]){
				print "error $val<=>$levels[$j-1] doesnot match $val<=>$hierarchy_check{$val} replacing with";
				my $tochange="$val;$levels[$j-1]";
				if(exists $hierarchy_change{$tochange}){
					$val=$hierarchy_change{$tochange};
					
				}else{
					$hierarchy_change{$tochange}= "$val $autonumber";
					$val=$hierarchy_change{$tochange};
					$autonumber++;
				}
				print " $val--(Unk)\n";
			}else{
				$hierarchy_check{$val} = $levels[$j-1];
			}
			$counts_line=$counts_line."$val\t";
		}else{
			
			if($j>=2 && exists $hierarchy_check{$val} && $hierarchy_check{$val} ne $levels[$j-1]){
				print "error $val<=>$levels[$j-1] doesnot match $val<=>$hierarchy_check{$val} replacing with";
				my $tochange="$val;$levels[$j-1]";
				if(exists $hierarchy_change{$tochange}){
					$val=$hierarchy_change{$tochange};
					
				}else{
					$hierarchy_change{$tochange}= "$val $autonumber";
					$val=$hierarchy_change{$tochange};
					$autonumber++;
				}
				print " $val--(else)\n";
			}else{
				$hierarchy_check{$val} = $levels[$j-1];
			}
			$counts_line=$counts_line."$val\t";
		}
	}
	shift(@entry);
	foreach(@entry){
		$counts_line=$counts_line."$_\t";
	}
	chop($counts_line);
	if($counts_line!~/^\#/){
	print $spf "$counts_line\n";
}
}
