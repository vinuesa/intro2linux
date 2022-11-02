#!/usr/bin/awk -f

#: PROGNAME: remove_zero_colsums_from_tsv.awk
#: AUTHOR: Pablo Vinuesa, CCG-UNAM; @pvinmex
#: FIRST VERSION: Nov 1st, 2022
#: AIM: remove columns from input table that have a zero colsum, for example abricate SUMMARY-resfinder.tsv,
#:       with [optional filtering providing a regex as first argument]. See EXAMPLES in Usage_Exit()

# Initialization and command-line parsing
BEGIN { 
   progname = "remove_zero_colsums_from_tsv.awk"
   version  = "0.2_2Nov22"  # v0.2_2Nov22 added the possibility to pass an [optional REGEX] as first argument
   
   # input is expected to be in tsv format; output also printed as tsv
   FS=OFS="\t"
      
   if (ARGC < 2) Usage_Exit(progname, version)

   if (ARGC == 3) # $0; regex; tsv file
   {
      # save the filtering regex in variable regex      
      regex = ARGV[1]
 
      # delete this first argument to avoid that,
      # in the main block below, awk treats it as a file 
      delete ARGV [1]
   }
}

# ACTION: Build datastuctures, filtering input line with regex, if provided as first argument on the command line
{ 
   for(i=1;i<=NF;i++) { 
      
      # parsing of input table if regex is provided as 1st argument
      if(regex !="" && NR == 1) {
          # 1. store complete line/record in the line AoA, field by field.
          line[NR][i]=$i
      }
      else if (regex !="" && $0 ~ regex && NR > 1 ){
          line[NR][i]=$i
     
          # 2. store column sums in col array
          col[i]+=$i
      }
      else if (regex !="" && $0 !~ regex && NR > 1 ) {
         # skip these non-matching lines from the data structures
	 next
      }
      
      # if no regex provided as 1st argument
      if (! regex){
          line[NR][i]=$i
     
          # 2. store column sums in col array
          col[i]+=$i
      }
   } 
}

# Print output line, if column count != 0
END {
 
   # 3. loop over lines/records; then over each record's column
   for ( l=1 ; l<=NR ; l++ )
   {
       printf line[l][1] "\t"

       # 4. loop over columns and print column if column sum != 0 (col[c] is defined)
       for (c=2;c<=NF;c++) if (col[c]) printf line[l][c] "\t"

       printf "\n"
    }
}

# Usage_Exit funtion
function  Usage_Exit(prog, vers) { 
   print "# USAGE FOR: " prog " v"vers > "/dev/stderr"
   print  "\t" prog " [<filtering string or regex (without //)>] <intput table (tsv) to filter>" > "/dev/stderr"
   print  "\n# AIM: Remove columns from input table that have a zero colsum" > "/dev/stderr"
   print  "\n# EXAMPLES:" > "/dev/stderr"
   print  "\t" prog " '(SBO|SPE)' SUMMARY-plasmidfinder.tsv | column -t" > "/dev/stderr"
   print  "\t" prog " SUMMARY-resfinder_SBO-SPE.tsv | column -t" > "/dev/stderr"
   exit
}

