blastparse.pl parse the m8 format into taxonomic profile of hits.

blast_parse_amr.pl
first will blast reads against a custom db (AMR genes in this case) and then create a counts file for each sample can be followed by parse_freq.pl to create a proper frequency profile of the present genes

retreive2.pl and dbmaker.pl
create a custom blast db using a entrez query. retreives sequences from NCBI automatically. requires:
fileconvert.pl
fasta_make.pl 
gb_utils_xml.pl

nanopore_dem_afterepi2me.pl
after runnin WIMP epi2me has stopped segregating reads into barcode folders. this script does that for us.

nanopore_epi2me_csv_parse scripts.
uses csv files created by epi2me for various functions like demultiplexing, Antimicrobial resistance, and also can use WIMP resuslts to create biom files.

virfindercommands.txt
referece commands for internal use.
