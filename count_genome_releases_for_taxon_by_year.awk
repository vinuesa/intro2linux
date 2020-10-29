# AUTHOR: Pablo Vinuesa, @pvinmex, https://www.ccg.unam.mx/~vinuesa/ @pvinmex
# source: https://github.com/vinuesa/intro2linux
# VERSION:0.1_2020-10-28

#>>> Initialize
BEGIN{ 
    # ARGC counts the arguments; ARGV[0] is awk
    if(ARGC < 2)     
    Usage_Exit()

    
    # we need to set FS to "\t" for propper parsing of assembly_summary.txt.gz, 
    # which is a tsv file
    FS="\t"
    print "taxon\trel_y\tn_genomes"
}

#>>> MAIN

# 1. filter taxon fields ($8) that match tax and susbstitute 2019/10/04 for 2019
$8 ~ tax {gsub(/\/.*$/, "")

# here we use a hash named count_y, using the year ($15) as key,
# autoimcrementing the associated values as we encounter them in each new record
count_y[$15]++}

# After filling the hash, print it out with some basic formatting
END{
      for (y in count_y) 
         if (y > 0) 
	   printf "%s\t%d\t%d\n", tax, y, count_y[y]
}


#>>> Function definition
function Usage_Exit() {
      # we call the function if the expected argument is not passed
      print "# USAGE:  zcat assembly_summary.txt.gz | awk -f count_genome_releases_for_taxon_by_year.awk tax=Pseudomonas";
      print "#   note1: need to pass tax='TAXON' as single argument at the end of the awk call as tax=<'TAXON STRING'>";
      print "#   note2: if using the gzip-compressed source file, you need to pipe it into awk with zcan, as shown above\n";
      print "# AIM: print the number of genomes released per year for a taxon provided as single argument,";
      print "#      based on data stored in assembly_summary.txt.gz";
      exit;
}

