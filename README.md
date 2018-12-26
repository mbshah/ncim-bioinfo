##Description of a few scripts in  this folder:##
---------------------------------------------------

KO2Path and Path2class:
https://github.com/mbshah/metgenomics

end-to-end-qiime.pl and end-to-end-qiime-swarm.pl
simple pipeline to carryout baisic qiime steps by just taking a few parametes as input.

merge.pl
will automatically merge mate pairs or paired end reads and also trim and carryout QC using trim gallore

present_or_not.pl
can convert profile tables to a table with only 0 or 1 to indicate either presence or absence.

split_sra_pairs.pl
split interleaved pairedends fasta files typically acquired from SRA into R1 and R2 files

runningt4f.pl
wrapper script for running tax4fun, required by end-to-end-qiime scipts.


##In virCodes folder:##
------------------------------
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
