#!/usr/bin/awk -f

#: PROGRAM: extract_CDSs_from_genbank.awk; first commmit Dec 23, 2020
#: AUTHOR: Pablo Vinuesa, @pvinmex, https://www.ccg.unam.mx/~vinuesa/
#: SOURCE: https://github.com/vinuesa/intro2linux
#: USAGE: extract_CDSs_from_genbank.awk genbank_file.gbk   
#: AIMS: 1. extracts the CDSs from a single GenBank file and writes them as input_genbank_basename_CDSs.fna
#:	 2. writes the GenBank file in tabular format as input_genbank_basename.tsv
#:	 3. extracts the complete DNA string from the input GenBank saving to input_genbank_basename.fsa
#: NOTE: the program assumes that the input GenBank file contains a single sequence record or LOCUS

BEGIN {   
    # Initializations
    DEBUG = 0  # set to 1 if debugging messages should be activated (prints to "/dev/stderr")
    start=end=complflag=pseudoflag=line_no=0
    product=locus_tag=""
    progname="extract_CDSs_from_genbank.awk"
    VERSION=0.2 # Dec 24, 2020; now captures full product name, even when split over two lines
                # Dec 23, 2020; first commit

    # check user input: needs the GenBank file to process
    #        as single argument on the command line
    if(ARGC < 2) 
        print_help(progname, VERSION)    
	
    # capture the gbk basename for naming output files written by this script
    input_gbk = ARGV[1]
    basename = input_gbk
    
    gsub(/\..*$/, "", basename)   
    gbk_tsv = basename".tsv"
    gbk_fsa = basename".fsa"
    gbk_CDSs = basename"_CDSs.fna"
    
    # hash used by rev_compl()
    compl_nucl["T"]="A"
    compl_nucl["A"]="T"
    compl_nucl["C"]="G"
    compl_nucl["G"]="C"
    compl_nucl["N"]="N"
}
#-------------------------------- END OF BEGIN BLOCK ----------------------------------#

# capture the locus ID on the first line to add to the FASTA header
# Then skip all lines, from the first record (line) to the source FEATURE,
# Note the use of an expr-REGEXP range: expr, /regexp/
NR == 1 { 
   locus_id = $2; 
   header   = locus_id
}

# skip lines until the source attribute of the FEATURES block is reached
NR == 1, /^\s{4,}source\s{4,}/ { next }

{
  # exclude pseudogenes, indicated by truncated gene coordinates with [<>]
  if (/^\s{4,}CDS\s{4,}[[:digit:]]+/ && flag == 0 && !/[<>]/ && !/complem/)
  {
     flag = 1 
     coord=$0   
     sub(/^\s{4,}CDS\s{4,}/, "", coord)
     #print "coord:"coord
     start=coord
     end=coord
     coord=0
     sub(/\.\.[[:digit:]]+$/, "", start)
     sub(/^[[:digit:]]+\.\./, "", end)
     
     complflag = 0
  }
  else if (/^\s{4,}CDS\s{4,}complement/ && ! flag && ! /[<>]/)
  {
     flag = 1
     coord=$0     
     sub(/^\s{4,}CDS\s{4,}complement\(/, "", coord)
     gsub(/[\)]/, "", coord)
     start = coord
     end = coord
     sub(/\.\.[[:digit:]]+/, "", start)
     sub(/^[[:digit:]]+\.\./, "", end)

     complflag = 1
  }
  
  if( flag && /^\s+\/locus_tag=/ )
  {
     locus_tag=$0
     sub(/^\s{8,}\/locus_tag=/, "", locus_tag)
     gsub(/["]/, "", locus_tag)
  } 

  # mark pseudogenes
  if( flag && /^\s{8,}\/pseudo/ ) pseudoflag = 1

  if(flag && /^\s{8,}\/product=/)
  {
     line_no = FNR
     product = $0
     sub(/^\s{8,}\/product=/, "", product)
     gsub(/["]/, "", product)

     if (pseudoflag) { product = product " [pseudogene]" }
     
     # we use geneID to generate a numeric index for the hash, so that it prints out in order!
     geneID++
     
     #printf "%s\t%s\t%d\t%d\t%d\n", locus_tag, product, start, end, complflag,  pseudoflag
     CDSs_string[geneID] = locus_tag "\t" product "\t" start "\t" end "\t" complflag "\t" pseudoflag
     CDSs_AoA[geneID]["l"]  = locus_tag
     CDSs_AoA[geneID]["p"]  = product
     CDSs_AoA[geneID]["s"]  = start
     CDSs_AoA[geneID]["e"]  = end
     CDSs_AoA[geneID]["c"]  = complflag
     CDSs_AoA[geneID]["pseudo"] = pseudoflag

     flag=start=end=complflag=revflag=pseudoflag=0
     locus_tag=""
  }
  
  # if the product description is split in two lines, 
  #     capture it and append it to the first line already saved in product
  if(FNR == line_no+1)
  {
       # this regexp is critical to ensure the caputre of
       #   all possible product names, like ANT(3'')-Ia, that must end with "
       # Set DEBUG = 1 at the beginning of the BEGIN{} block
       #   to print all the second lines for producut names to "/dev/stderr"
       #if($0 ~/^\s+[[:alnum:] ]+\"$/)
       if($0 ~/^\s+[a-zA-Z0-9\(\)\'\- ]+"$/) # <== matches patterns like 'gene ANT(3'')-Ia'
       {
	   prod_2cnd_line = $0
	   gsub(/^\s+/, "", prod_2cnd_line)
	   gsub(/["]/, "", prod_2cnd_line)
	   product = product " " prod_2cnd_line
	   
	   if(DEBUG)
	      print "FNR=" FNR, "line_no="line_no, "prod_2cnd_line="prod_2cnd_line > "/dev/stderr"

	   CDSs_AoA[geneID]["p"] = product
	   line_no=product=""
       }
  }
  
  # After the ORIGING mark, starts the dna string
  if (flag == 0 && !/^ORIGIN/) { 
      next 
  }
   
  # Remove spaces and digits preceding sequence strings
  #  remove also the single spaces between 10 bp sequence strings
  if (/^ORIGIN/) { flag = 1; next}
  
  if(flag && $0 ~ /^\s+[[:digit:]]+\s+[AGCTNagctn\s]+/)
  {   
     sub(/^\s+[[:digit:]]+\s+/, "", $0)
     gsub(/\/\//, "", $0)
     gsub(/[[:space:]]/, "", $0)

     seq= seq $0
  }
}
# ------------------------------- END OF PATTERN & ACTION BLOCK -----------------------------------#

END {
    # save the genbank's DNA string in the dna hash, indexed by LOCUS
    seq = toupper(seq)
    dna[header] = seq
    	
    # conditional required to avoid errors printed to STDERR when the help menu is printed
    if( ARGC > 1 ) 
    {  
    	# Need to remove pre-existing "basename".tsv" file, as the script appends lines to it 
    	#     (print stuff >> gbk_tsv)
        # Can use any of the following options:
    	#    1. call system()
    	#    2. pipe command to sh
    	#    the second may by more portable and is used as a function call
        #
    	#if(system("[ -s " gbk_tsv " ]") == 0) { system("rm "gbk_tsv) }
    	#OR:
    	#rm_if_exists="[ -s " gbk_tsv " ] && rm "gbk_tsv
    	#print rm_if_exists | "/bin/sh"
    	#close("/bin/sh")
        rm_if_exists(gbk_tsv)
    
    	# print GBK table header for gbk_tsv
    	print "CDS_no\tlocus_tag\tproduct_description\tstart\tend\tcomplement\tpseudogene" > gbk_tsv
    
    
    	# if exists, rm the FASTA file holding CDSS, since the script appends to it
    	#     (print stuff >> gbk_CDSs)
    	rm_if_exists(gbk_CDSs)
    
    	# loop over the CDSs_AoA array of arrays to print 
    	#      i) the GBK file in tabular format to file gbk_tsv
    	#      ii) and the CDSs to STDOUT
    	for (k in CDSs_AoA)
    	{
            # 1. print the GBK file in tabular format to file gbk_tsv
            printf "%d\t%s\t%s\t%d\t%d\t%d\t%d\n", k, CDSs_AoA[k]["l"], CDSs_AoA[k]["p"], CDSs_AoA[k]["s"], 
        	  CDSs_AoA[k]["e"], CDSs_AoA[k]["c"], CDSs_AoA[k]["pseudo"] >> gbk_tsv
    	
            # 2. print the CDSs to STDOUT
            extract_sequence_by_coordinates(CDSs_AoA[k]["l"]" "CDSs_AoA[k]["p"], seq, CDSs_AoA[k]["s"], 
        	  CDSs_AoA[k]["e"], CDSs_AoA[k]["c"], gbk_CDSs) 
    	}
        # report if the files were successfully written to disk
    	print_if_exists(gbk_tsv)
    	print_if_exists(gbk_CDSs)
    	 
    	# if exists, rm the FASTA file holding the DNA string, since the script appends to it
    	#     (print stuff >> gbk_fsa)
        rm_if_exists(gbk_fsa)

    	# print the complete DNA sequence string of the input GenBank to to file gbk_fsa
    	for (h in dna)
    	{
             printf ">%s\n%s\n", header, seq >> gbk_fsa
    	}
        # report if the file was successfully written to disk
        print_if_exists(gbk_fsa)
    }
}
# --------------------------------------------- END OF END BLOCK ---------------------------------------- #

#=========================================================================================================#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION DEFINITIONS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
#---------------------------------------------------------------------------------------------------------#

function rm_if_exists(file,       cmd)
{
    cmd="[ -s " file " ] && rm "file
    print cmd | "/bin/sh"
    close("/bin/sh")
}
#---------------------------------------------------------------------------------------------------------

function print_if_exists(file,     cmd)
{
    cmd="[ -s " file " ] && echo '# wrote file '"file
    print cmd | "/bin/sh"
    close("/bin/sh")
}
#---------------------------------------------------------------------------------------------------------

function rev_string(seq       ,i,s)
{
   for(i=length(seq); i != 0; i--)
   {
      s = s substr($0, i, 1)
   }
   return s
}
#---------------------------------------------------------------------------------------------------------

function rev_compl(header, dna_seq, outfile,       i, k, s) # reverse complement
{  # receives two arguments: the fasta header and the CDS sequence to be translated
    i=k=s=""
    
    dna_seq = toupper(dna_seq) # to match the keys of compl_nucl
    for( i = length(dna_seq); i !=0; i-- ) { # note that the loop reads the sequence from end to beginnig
         k = substr(dna_seq, i, 1)
         s = s compl_nucl[k]
    }
    printf ">%s\n%s\n", header, s >> outfile
}
#---------------------------------------------------------------------------------------------------------

function extract_sequence_by_coordinates(header, dna_seq, start, end, revCompFlag, outfile,        seq)
{   
    if(start == 1 && revCompFlag)
    {
    	seq = substr(dna_seq, start, end)
	rev_compl(header, seq, outfile)
    }

    if(start == 1 && ! revCompFlag)
    {
    	print ">"header"\n"substr(dna_seq, start, end) >> outfile
    }		 
    
    if( start > 1 && ! revCompFlag ) 
    {
        print ">"header"\n"substr(dna_seq, start, end-start+1) >> outfile
    }
    
    if( start > 1 && revCompFlag ) 
    {
	seq = substr(dna_seq, start, end-start+1)
	rev_compl(header, seq, outfile)
    }
}
#---------------------------------------------------------------------------------------------------------

function print_help(prog, vers) 
{
  print "# AIMS:" > "/dev/stderr" 
  print "#  1. extracts the CDSs from a single GenBank file, saving them as input_genbank_basename_CDSs.fna" > "/dev/stderr"
  print "#  2. writes the GenBank file in tabular format as input_genbank_basename.tsv" > "/dev/stderr"
  print "#  3. extracts the complete DNA string from the input GenBank, saving it as input_genbank_basename.fsa" > "/dev/stderr"
  print "# NOTE", prog, "v"vers, "assumes that the input GenBank file contains a single sequence record or LOCUS" > "/dev/stderr"
  print "# USAGE of", prog, "v"vers":" > "/dev/stderr"
  printf "\t%s %s\n\n", prog, "genbank_file_with_single_record.gbk" > "/dev/stderr"  
  
  exit 0
}
