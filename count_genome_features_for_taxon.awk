# AUTHOR: Pablo Vinuesa, CCG-UNAM; https://www.ccg.unam.mx/~vinuesa/; twitter: @pvinmex
# source: https://github.com/vinuesa/intro2linux
# VERSION:0.1_2020-11-1
# AIM: get summary statistics (counts) of the features found in user-provided field and taxon
#      in NCBI's assembly_summary.txt.gz genome assembly table
# ---------------------- #
# >>> Initialization <<< #
# ---------------------- #
BEGIN { 
    # ARGC counts the arguments; ARGV[0] is awk
    if(ARGC < 3)   # needs two positional arguments 
       Usage_Exit()
   
    # we need to set FS"\t" for propper parsing of assembly_summary.txt.gz, 
    # which is a tsv file; set also OFS=FS
    FS=OFS="\t"
        
    # Check that user provides correct arguments
    #   use sub() to substitute the 'tax=' and 'f=' prefixes for ""
    #   to retain only the option values passed on the command line   
    tax = ARGV[1]
    if( sub(/tax=/, "", tax) ) {       
        sub(/tax=/, "", tax)
    } 
    else {
       print "ERROR: you need to provide tax=TAXON as the first argument to the script!"
       print "       ==> you provided", tax
       print ""
       Usage_Exit()   
    }
    
    idx = ARGV[2]
    if ( sub(/f=/, "", idx) ) {
         sub(/f=/, "", idx)
    }
    else {
       print "ERROR: you need to provide f=<table_field_number> as the second argument to the script!"
       print "       ==> you provided", idx
       print ""
       Usage_Exit()   
    }
    
    # check that meaningful column indexes are passed to the script
    if ( idx !~ /^(2|5|11|12|13|14|15)$/ ) { Usage_Exit() }
        
    # sort array numeric output values in ascending order
    PROCINFO["sorted_in"] = "@val_num_asc"
}

# ------------ #
# >>> MAIN <<< #
# ------------ #
# 1. fill the hash features with the field names, indexed by consecutive integers/positions
#      using the split() function applied to record # 2 (the column header) 
#      to add the header field names to consecutive numeric indices such as:
#      features[1] = "assembly_accession"; features[2] = "bioproject"; ...
NR==2 { split( $0, features, "\t", seps ) } 

# 2. filter taxon fields ($8) that match /tax/ 
$8 ~ tax { 
    # globally susbstitute months and days by nothing; 
    # i.e: change 2019/10/04 for 2019
    gsub( /\/.*$/, "" )

    # we use a hash named count_f, using the feature index ($idx) as the hash's key,
    # to hold counts of the associated values as we encounter them in each new record
    count_f[$idx]++
}

# 3. END{} block. print the hash values in count_f indexed by f
END {
       #for (i in features)  { print i, features[i] }
       if ( idx ~ /^(2|5|11|12|13|14|15)$/ ) {
           print "taxon", features[idx], "n_genomes"
           for ( f in count_f ) 
              if ( f > 0 || f != "" ) 
	         print tax, f, count_f[f]
       }	         
}

# --------------------------- #
# >>> Function definition <<< #
# --------------------------- #
function Usage_Exit() {
      # the function is called if the proper arguments are not passed to the script
      print "# USAGE: zcat assembly_summary.txt.gz | awk -f count_genome_features_for_taxon.awk tax=Pseudomonas f=12"
      print "#    note1: need to pass tax='TAXON' and f=<'field number'[2|5|11-15]> as positional arguments at the end of the awk call, as shown" 
      print "#    note2: if using the gzip-compressed source file, you need to pipe it into awk with zcat, as shown above" 
      print "# AIM: print a table with the the number of features found in assembly_summary.txt.gz for a taxon" 
      print "#      both provided as arguments to the program in the format tax=TAXON f=feature_field_number" 
      exit;
}

