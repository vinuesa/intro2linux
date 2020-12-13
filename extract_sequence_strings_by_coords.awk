# extract_sequence_strings_by_coords.awk
# AUTHOR: Pablo Vinuesa, @pvinmex, https://www.ccg.unam.mx/~vinuesa/
# source: https://github.com/vinuesa/intro2linux
# USAGE:
#   awk -f extract_sequence_strings_by_coords.awk start=2 end=5 seq.fasta

BEGIN {
     progname = "extract_sequence_strings_by_coords.awk"
     version  = 0.1  # dec 12, 2020
    
     if  (ARGC < 4) Usage_Exit(progname, version)
    
     RS=">"
     FS="\n"
     
     # check that user provides proper arguments
     if(ARGV[1] !~ /^start=/) Usage_Exit(progname, version)
     if(ARGV[2] !~ /^end=/) Usage_Exit(progname, version)
     if(ARGC == 5 && ARGV[3] !~ /^name=/) Usage_Exit(progname, version)  
}

NR > 1 {  
   seq="" 
   for (i=2;i<=NF;i++) seq=seq""$i
       
   if(ARGC == 5)
   {
      #print "name=", name, "$1=", $1
      if ( $1 == name && start == 1) print ">"$1"\n"substr(seq,start,end)
      if ( $1 == name && start > 1 ) print ">"$1"\n"substr(seq,start,end-start+1)
   }
   else
   {
       if(start == 1) print ">"$1"\n"substr(seq,start,end)
       if(start > 1)  print ">"$1"\n"substr(seq,start,end-start+1)
   }
}  

# function definition
function Usage_Exit(prog, vers) {
 
   print "# USAGE FOR", prog, "v"vers > "/dev/stderr"
   print "#  ", prog, "start=<int> end=<int> [name=seq_name] fasta_file" > "/dev/stderr"
   print "# Example: awk -f", prog, "start=2 end=5 name=seq1 seq.fasta" > "/dev/stderr"
   print "# NOTES: 1. the fasta_file can contain 1 or more sequences" > "/dev/stderr"
   print "#        2. the third argument name=sequence_name is optional." > "/dev/stderr"
   print "#            Use it to filter out a specific sequence from a multifasta" > "/dev/stderr"
   exit;
}
