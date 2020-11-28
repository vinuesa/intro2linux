#!/usr/bin/gawk -f

#: AUTHOR: Pablo Vinuesa, @pvinmex, https://www.ccg.unam.mx/~vinuesa/
#: source: https://github.com/vinuesa/intro2linux
#: USAGE: extract_DNA_string_from_genbank.awk genbank_file.gbk   
#: AIM: extracts the DNA string from a GenBank file and prints to STDOUT
#: NOTE: the program assumes that the input GenBank file contains a single sequence record or LOCUS

BEGIN {

   progname="extract_DNA_string_from_genbank.awk"
   VERSION=0.1  # Nov 27, 2020

   if(ARGC < 2) # needs the genbank file to process as single argument on the command line
        Usage_Exit()    
	
   # print the gbk basename as the FASTA header
   input_gbk=ARGV[1]
   gsub(/\..*$/, "", input_gbk)   
   print ">"input_gbk
}

# skip all lines, from the first record (line) to the ORIGIN flag,
#   which precede the sequnce string to be extracted
# Note the use of an expr-REGEXP range: expr, /regexp/
NR == 1, /^ORIGIN/ { next } 

# remove spaces and digits preceding sequence strings
#  remove also the single spaces between 10 bp sequence strings
{ sub(/^\s+[[:digit:]]+\s+/, ""); gsub(/\/\//, ""); gsub(/[[:space:]]/, "") }

# skip any empty lines and convert sequence strings to uppercase
NF > 0 { print toupper($0)}

function Usage_Exit() {
  
  print "# AIM: extracts the DNA string from a GenBank file and prints to STDOUT"
  print "# NOTE", progname, "v"VERSION, "assumes that the input GenBank file contains a single sequence record or LOCUS"  
  print "# USAGE of", progname, "v"VERSION":" 
  printf "\t%s %s\n\n",  progname, "genbank_file.gbk"   
  
  exit 0

}
