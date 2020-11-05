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
         version  = 0.2  # nov 04, 2020
         
         if  (ARGC < 2) Usage_Exit(progname, version)

         # Capturamos la cadena que vamos a usar para filtrar
         s = ARGV [1]
         
         # y borramos el primer argumento para que awk no lo trate como un archivo
         #  e intente abrirlo. 
         delete ARGV [1];
         
	 RS=">"
}

# MAIN; if line matches filtering string, print record ;)
$0~s{print ">"$0}

# function definition
function  Usage_Exit  (prog, vers) {
 
   print "# USAGE FOR", prog, "v"vers > "/dev/stderr"
   print  p, "<filtering_string>  <multifasta_file>" > "/dev/stderr"
   print "#   Pass a string as first argument to filter the FASTA_file," > "/dev/stderr"
   print "#   provided as second argument, printing only records matching the string" > "/dev/stderr"
   exit;
}
