#!/usr/bin/awk -f

# AUTHOR: Pablo Vinuesa, @pvinmex, https://www.ccg.unam.mx/~vinuesa/
# source: https://github.com/vinuesa/intro2linux
# VERSION:0.1
# Usage:  filter_fasta_sequences.awk  <filtering_string>  <multifasta_file>
#   Read a string from STDIN to filter fasta_file, 
#   to print out only the subset of the fasta_file 
#   matching the filtering string


BEGIN   {
         if  (ARGC < 1)
             Usage_Exit();


         s = ARGV [1];
         delete ARGV [1];


	 RS=">"
        }


# MAIN; if line matches filtering string, print record ;)
$0~s{print ">"$0}


# function definition
function  Usage_Exit  ()
  {
   print "# Usage:  filter_fasta_sequences.awk <filtering_string>  <multifasta_file>";
   print "#   Read a string from STDIN to filter fasta_file,";
   print "#   to print out only the subset of the fasta_file";
   exit;
  }
