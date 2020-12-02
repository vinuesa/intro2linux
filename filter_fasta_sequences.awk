#!/usr/bin/awk -f

# AUTHOR: Pablo Vinuesa, @pvinmex, https://www.ccg.unam.mx/~vinuesa/
# source: https://github.com/vinuesa/intro2linux
# Usage:  filter_fasta_sequences.awk  <filtering_string>  <multifasta_file>
#   Read a string from STDIN to filter fasta_file, 
#   to print out only the subset of the fasta_file 
#   matching the filtering string
# NOTE: this is a demo script to teach basic awk programming

BEGIN {
         progname = "filter_fasta_sequences.awk"
         version  = 0.3  # dec 02, 2020
         
         if  (ARGC < 2) Usage_Exit(progname, version)

         # save te filtering string in variable s
         s = ARGV [1]
         
         # delete this first argument with delete.
	 # This is to avoid that in the main block below, 
	 # the command interpreter treats it as a file 
         delete ARGV [1];
         
	 RS=">"
}

# MAIN; if line matches filtering string, print record ;)
$0 ~ s { print ">"$0 }

# function definition
function  Usage_Exit(prog, vers) {
 
   print "# USAGE FOR", prog, "v"vers > "/dev/stderr"
   print  prog, "<filtering_string>  <multifasta_file>" > "/dev/stderr"
   print "#   Pass a string as first argument to filter the FASTA_file," > "/dev/stderr"
   print "#   provided as second argument, printing only records matching the string" > "/dev/stderr"
   exit;
}
