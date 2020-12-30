#!/usr/bin/awk -f

#---------------------------------------------------------------------------------------------------------
#: PROGRAM: extract_CDSs_from_GenBank.awk; first commmit Dec 23, 2020
#: AUTHOR: Pablo Vinuesa, @pvinmex, https://www.ccg.unam.mx/~vinuesa/
#: SOURCE: https://github.com/vinuesa/intro2linux
#: USAGE: extract_CDSs_from_genbank.awk genbank_file.gbk   
#: AIMS: 
#:   1. extracts the CDSs from a single GenBank file and writes them to genbank_basename_CDSs.fna
#:   2. translates the CDSs and writes them to file genbank_basename_proteome.faa
#:   3. writes the GenBank file in tabular format as genbank_basename.tsv
#:   4. extracts the complete DNA string(s) from the input GenBank saving them to genbank_basename.fsa
#: NOTES: 
#:   1. the program currently does not deal with CDSs containing introns.
#:       CDSs containing the join statement: complement(join(4497616..4498557,4498557..4498814)), are skipped
#:       Use only for bacterial|mitochondrial|plastid|yeast genomes
#=========================================================================================================

BEGIN {   
    # Initializations
    DEBUG = 1  # set to 1 if debugging messages should be activated (prints to "/dev/stderr")
    
    # print the FASTA sequences stored in hashes in ascending order by locus_tag number 
    PROCINFO["sorted_in"] = "@ind_num_asc"
    
    progname="extract_CDSs_from_GenBank.awk"
    VERSION=0.6 # Dec 29, 2020; 
                #  * now also captures the protein_id to label prots in *faa files; it is added to tsv outfile
                #  * Only CDSs are labeled with locus_tag
                #  * Manages input gbks like E. coli MT559985 that do not have locus_tag attributes,
		#      labeling them with locus_id"_CDS"geneID
       # v0.5 Dec 28, 2020; slightly more streamlined (awkish) code, by using more awk defaults (avoid $0 ~ /regexp/, etc)
       #    prints LOCUS_length to *tsv file using a nested AoA as the 1ary dna AoA key: dna[CDSs_AoA[k]["LOCUS"]]["len"]
       # v0.4 Dec 26, 2020; prints also the CDS's translation products to file
       # v0.3 Dec 25, 2020; correctly appends the [pseudogene] label to the end of the product line
       # v0.2 Dec 24, 2020; now captures full product name, even when split over two lines
       # v0.1 Dec 23, 2020; first commit

    # check user input: needs the GenBank file to process
    #        as single argument on the command line
    if(ARGC < 2) 
        print_help(progname, VERSION)    
	
    # capture the gbk basename for naming output files written by this script
    input_gbk = ARGV[1]
    basename = input_gbk
    
    gsub(/\..*$/, "", basename)   
    gbk_tsv = basename".tsv"
    gbk_faa = basename"_proteome.faa"
    gbk_fsa = basename".fsa"
    gbk_CDSs = basename"_CDSs.fna"
    
    # hash used by rev_compl()
    compl_nucl["T"]="A"
    compl_nucl["A"]="T"
    compl_nucl["C"]="G"
    compl_nucl["G"]="C"
    compl_nucl["N"]="N"

    # initialize the "codons" hash holding the codon-aminoacid pairs, 
    #	based on the universal genetic code, to translate CDSs with translate_dna()
    codons["ATA"]="I"; codons["ATC"]="I"; codons["ATT"]="I"; codons["ATG"]="M";
    codons["ACA"]="T"; codons["ACC"]="T"; codons["ACG"]="T"; codons["ACT"]="T";
    codons["AAC"]="N"; codons["AAT"]="N"; codons["AAA"]="K"; codons["AAG"]="K";
    codons["AGC"]="S"; codons["AGT"]="S"; codons["AGA"]="R"; codons["AGG"]="R";
    codons["CTA"]="L"; codons["CTC"]="L"; codons["CTG"]="L"; codons["CTT"]="L";
    codons["CCA"]="P"; codons["CCC"]="P"; codons["CCG"]="P"; codons["CCT"]="P";
    codons["CAC"]="H"; codons["CAT"]="H"; codons["CAA"]="Q"; codons["CAG"]="Q";
    codons["CGA"]="R"; codons["CGC"]="R"; codons["CGG"]="R"; codons["CGT"]="R";
    codons["GTA"]="V"; codons["GTC"]="V"; codons["GTG"]="V"; codons["GTT"]="V";
    codons["GCA"]="A"; codons["GCC"]="A"; codons["GCG"]="A"; codons["GCT"]="A";
    codons["GAC"]="D"; codons["GAT"]="D"; codons["GAA"]="E"; codons["GAG"]="E";
    codons["GGA"]="G"; codons["GGC"]="G"; codons["GGG"]="G"; codons["GGT"]="G";
    codons["TCA"]="S"; codons["TCC"]="S"; codons["TCG"]="S"; codons["TCT"]="S";
    codons["TTC"]="F"; codons["TTT"]="F"; codons["TTA"]="L"; codons["TTG"]="L";
    codons["TAC"]="Y"; codons["TAT"]="Y"; codons["TAA"]="*"; codons["TAG"]="*";
    codons["TGC"]="C"; codons["TGT"]="C"; codons["TGA"]="*"; codons["TGG"]="W";
    
    unknown = "X"
}
#-------------------------------- END OF BEGIN BLOCK ----------------------------------#

# The following block will parse the relevant annotation features of CDSs,
#   saving them in the the CDSs_AoA (Array of Arrays)

# Capture the locus ID on the first line to add to the FASTA header
# Then skip all lines, from the first record (line) to the source FEATURE,
# Note the use of an expr-REGEXP range: expr, /regexp/
FNR == 1 || /^LOCUS/ { 
    locus_id = $2 
   
    printf(">>> Processing LOCUS %s of %s ...\n", locus_id, FILENAME) > "/dev/stderr"
}

# skip lines until the source attribute of the FEATURES block is reached
/^LOCUS/, /^\s{4,}source\s{4,}/ { next }

{
  # Exclude pseudogenes, indicated by truncated gene coordinates with [<>], 
  #  or by the join statement: complement(join(4497616..4498557,4498557..4498814))
  if (/^\s+CDS\s{4,}[[:digit:]]+/ && flag == 0 && !/[<>,]/ && !/complem/ && !/join/)
  {
     flag  = 1 
    
     # note the pattern capturing parentheses used in this regexp, 
     #   which capture the coords 123..789 into the c_arr array
     match($0, /^\s{4,}CDS\s{4,}([[:digit:]]+)\.\.([[:digit:]]+)/, c_arr)
     start = c_arr[1]
     end   = c_arr[2]

     complflag = 0
  }
  
  if (/^\s+CDS\s{4,}complement/ && ! flag && ! /[<>,]/  && !/join/)
  {
     flag  = 1
     coord = $0     
     
     # for the sake of showing different string-manipulation functions
     #   in this case we combine sub(), gsub() and split()
     #   to capture the CDS's coords
     sub(/^\s{4,}CDS\s{4,}complement\(/, "", coord)
     gsub(/[\)]/, "", coord)

     split(coord, c_arr, /\.\./)
     start = c_arr[1]
     end   = c_arr[2]
     
     complflag = 1
  }
  
  if( flag && /^\s+\/locus_tag=/ )
  {
     match($0, /^\s{8,}\/locus_tag="([[:alnum:]_]+)"/, lt_arr)
     locus_tag = lt_arr[1]
  } 
    
  # mark pseudogenes
  if( flag && /^\s{8,}\/pseudo/ ) pseudoflag = 1

  if(/^\s{8,}\/product=/)
  {
     line_no = FNR
     product = $0
     sub(/^\s{8,}\/product=/, "", product)
     gsub(/["]/, "", product)
     gsub(/[>]/, "=", product)  # # remove '>' from (adenine(2058)-N(6))-dimethyltransferase > Erm(B)"
     
     # we use geneID to generate a numeric index for the hash, so that it prints out in order!
     geneID++
     
     # All the CDS's-relevant information parsed from the annotations 
     #   will be stored in the CDSs_AoA (Array of Arrays)
     if (locus_tag) CDSs_AoA[geneID]["l"] = locus_tag
     
     # some GenBank files like Escherichia_coli MT559985 do not have a locus_tag!
     if (!locus_tag) CDSs_AoA[geneID]["l"] = locus_id"_CDS"geneID
     CDSs_AoA[geneID]["LOCUS"]  = locus_id
     CDSs_AoA[geneID]["s"]      = start
     CDSs_AoA[geneID]["e"]      = end
     CDSs_AoA[geneID]["c"]      = complflag
     CDSs_AoA[geneID]["pseudo"] = pseudoflag

     complflag=0
  }

  # if the product description is split in two lines, 
  #     capture it and append it to the first line already saved in product
  if(FNR == line_no+1)
  {
     # this regexp is critical to ensure the caputre of
     #   all possible product names, like ANT(3'')-Ia, or
     #   (adenine(2058)-N(6))-dimethyltransferase > Erm(B)", that must end with "
     # Set DEBUG = 1 at the beginning of the BEGIN{} block
     #   to print all the second lines for producut names to "/dev/stderr"
     #if($0 ~/^\s+[[:alnum:] ]+\"$/)
     if( /^\s+[a-zA-Z0-9\(\)\'\-\[\],;:> ]+"$/ )
     {
     	prod_2cnd_line = $0
     	gsub(/^\s+/, "", prod_2cnd_line)
     	gsub(/["]/, "", prod_2cnd_line)
	gsub(/[>]/, "=", prod_2cnd_line) # (adenine(2058)-N(6))-dimethyltransferase > Erm(B)"
     	product = product " " prod_2cnd_line
    	
     	if(DEBUG)
     	   print "FNR=" FNR, "line_no="line_no, "locus_tag="locus_tag, 
	          "prod_2cnd_line="prod_2cnd_line > "/dev/stderr"

     	if (pseudoflag) product = product " [pseudogene]"

     	CDSs_AoA[geneID]["p"] = product
     	line_no=pseudoflag=0
     }
     else
     {
     	if (pseudoflag) product = product " [pseudogene]"

     	CDSs_AoA[geneID]["p"] = product
     	line_no=pseudoflag=0
     }
  }
  
  if( /^\s+\/protein_id="/ )
  {
     match($0, /^\s{8,}\/protein_id="(\S+)"$/, pi_arr)
     CDSs_AoA[geneID]["pid"] =  pi_arr[1]

     if (DEBUG) print "protein_id="pi_arr[1]
     flag = 0
  } 

  # After the ORIGING mark, starts the dna string
  if (flag == 0 && !/^ORIGIN/)  
      next 
   
  # Remove spaces and digits preceding sequence strings
  #  remove also the single spaces between 10 bp sequence strings
  if (/^ORIGIN/) { flag = 1; seq = ""; next }
  
  if(flag && /^\s+[[:digit:]]+\s+[AGCTNagctn\s]+/)
  {   
     sub(/^\s+[[:digit:]]+\s+/, "")
     gsub(/[[:space:]]/, "")

     seq = seq $0
  }

  # fill the dna AoA, once we have reached the end of record mark '//'
  if(flag && /^\/\/$/)
  {
     # save the genbank's DNA string and lenght in the dna array of arrays, 
     # indexed by [locus_id]["seq"] and [locus_id]["len"], respectively
     seq = toupper(seq)
     dna[locus_id]["seq"] = seq
     dna[locus_id]["len"] = length(seq)

     flag = 0
     locus_id=""
  }
}
# ------------------------------- END OF PATTERN-ACTION BLOCK -----------------------------------#

END {
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
       print "CDS_no\tLOCUS\tLOCUS_length\tlocus_tag\tprotein_id\tproduct_description\tstart\tend\tcomplement\tpseudogene" > gbk_tsv
    
       # if exists, rm the FASTA file holding CDSS, since the script appends to it
       #     (print stuff >> gbk_CDSs)
       rm_if_exists(gbk_CDSs)
    
       # loop over the CDSs_AoA array of arrays to print 
       #      i) the GBK file in tabular format to file gbk_tsv
       #      ii) and the CDSs to STDOUT
       for (k in CDSs_AoA)
       {
           # 1. print the GBK file in tabular format to file gbk_tsv; Nothe the use of an AoA as 1ary key to the dna AoA
           printf "%d\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n", k, CDSs_AoA[k]["LOCUS"], dna[CDSs_AoA[k]["LOCUS"]]["len"],
	            CDSs_AoA[k]["l"], CDSs_AoA[k]["pid"], CDSs_AoA[k]["p"], CDSs_AoA[k]["s"], CDSs_AoA[k]["e"], 
		      CDSs_AoA[k]["c"], CDSs_AoA[k]["pseudo"] >> gbk_tsv
    
           # 2. extract_sequence_by_coordinates, write CDSs.fna, translate CDSs and write proteome.faa
           extract_sequence_by_coordinates(CDSs_AoA[k]["l"], CDSs_AoA[k]["pid"], CDSs_AoA[k]["p"], seq, CDSs_AoA[k]["s"], 
        	 CDSs_AoA[k]["e"], CDSs_AoA[k]["c"], gbk_CDSs, CDSs_AoA[k]["pseudo"]) 
       }
       # report if the files were successfully written to disk
       print_if_exists(gbk_tsv)
       print_if_exists(gbk_CDSs)
        
       # if exists, rm the FASTA file holding the DNA string, since the script appends to it
       #     (print stuff >> gbk_fsa)
       rm_if_exists(gbk_fsa)

       # print the complete DNA sequence strings for each LOCUS/record 
       #  in the input GenBank to file gbk_fsa
       for (h in dna)
       {
            printf ">%s\n%s\n", h, dna[h]["seq"] >> gbk_fsa
       }
       # report if the file was successfully written to disk
       print_if_exists(gbk_fsa)

       # translate fastas and print to disk
       rm_if_exists(gbk_faa)

       for (h in prots) {
           # make sure we print only non-redundant sequences
           hcount[h]++
           
           # skip printing any pesudogene translations, which will be empty
           if( h ~ / \[pseudogene\]$/ ) continue
           
           printf(">%s\n%s\n", h, prots[h]) >> gbk_faa
       } 
       # report if the file was successfully written to disk
       print_if_exists(gbk_faa)
    }
}
# ------------------------------------------- END BLOCK FINISHES ---------------------------------------- #

#=========================================================================================================#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION DEFINITIONS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
#---------------------------------------------------------------------------------------------------------#

function rm_if_exists(file,       cmd)
{
    cmd="[ -s " file " ] && rm "file
    #print cmd | "/bin/sh" # <-- symbolic link: bin/sh -> dash
    #close("/bin/sh")
    print cmd | "/usr/bin/env bash"
    close("/usr/bin/env bash")
}
#---------------------------------------------------------------------------------------------------------

function print_if_exists(file,     cmd)
{
    cmd="[ -s " file " ] && echo '# wrote file '"file
    print cmd | "/usr/bin/env bash"
    close("/usr/bin/env bash")
}
#---------------------------------------------------------------------------------------------------------

function rev_string(seq,       i,s)
{
   for(i=length(seq); i != 0; i--)
   {
      s = s substr($0, i, 1)
   }
   return s
}
#---------------------------------------------------------------------------------------------------------

function rev_compl(lt, pi, pro, dna_seq, outfile,       i, k, s, dna_header, prot_header) # reverse complement
{  # receives two arguments: the fasta header and the CDS sequence to be translated
    i=k=s=""
    
    dna_header  = lt " " pro
    prot_header = pi " " pro

    dna_seq = toupper(dna_seq) # to match the keys of compl_nucl
    for( i = length(dna_seq); i !=0; i-- ) { # note that the loop reads the sequence from end to beginnig
         k = substr(dna_seq, i, 1)
         s = s compl_nucl[k]
    }
    printf(">%s\n%s\n", dna_header, s) >> outfile
    translate_dna(prot_header, s)
}
#---------------------------------------------------------------------------------------------------------

function extract_sequence_by_coordinates(ltg, pid, prod, dna_seq, start, end, revCompFlag, outfile,        seq, dna_header, prot_header)
{   
    # This function extracts the CDS string using the caputred coordinates
    #  and appends the sequence to the output FASTA file, calling the auxiliary
    #  function rev_compl() if required. 
    #  Both extract_sequence_by_coordinates() and rev_compl() call translate_dna()
    #    to translate the CDSs using the codons hash, saving the whole proteome in 
    #    the prots hash
    dna_seq = toupper(dna_seq)
    
    dna_header  = ltg " " prod
    prot_header = pid " " prod
    
    if(start == 1 && revCompFlag)
    {
    	seq = substr(dna_seq, start, end)
	rev_compl(ltg, pid, prod, seq, outfile)
    }

    if(start == 1 && ! revCompFlag)
    {
    	print ">"dna_header"\n"substr(dna_seq, start, end) >> outfile
	translate_dna(prot_header, substr(dna_seq, start, end))
    }		 
    
    if( start > 1 && ! revCompFlag ) 
    {
        print ">"dna_header"\n"substr(dna_seq, start, end-start+1) >> outfile
	translate_dna(prot_header, substr(dna_seq, start, end-start+1))
    }
    
    if( start > 1 && revCompFlag ) 
    {
	seq = substr(dna_seq, start, end-start+1)
	rev_compl(ltg, pid, prod, seq, outfile)
    }
}
#---------------------------------------------------------------------------------------------------------

function translate_dna(header, dna_seq,       s, i, p)
{  # receives two arguments: the fasta header and the CDS sequence to be translated
   
   # Initialize variables: 
   #  do-while loop control variable i (nt counter) 
   #   and p, which will hold the translation product
   {i=1; p=""; triplet_counter=0}

   # Here we run a do-while loop; the do loop is a variation of the while looping statement. 
   #  The do loop always executes the body once and then repeats the body as long as the condition is true
   # We use the do-while loop, to get a first triplet string saved in s; 
   #  then the while loop keeps going until substr() got the last triplet, resulting in an empty s="".
   do {
        # First check that the script got some input
        #   if not, exit with an error message
        if(length(dna_seq) == 0) {
            print "ERROR: need a DNA sequence string to translate (valid DNA sequence, divisible by 3)" > "/dev/stderr"
            exit 1
       
        # Check that the DNA sequence string is divisible by 3 using the modulo operator
        #   if not, exit with an error message
        } else if(length(dna_seq)%3) { 
            printf("# WARNING: input DNA sequence for %s not divisible by 3. Will skip it!\n", header) > "/dev/stderr"
            continue
        }
        
        # use substr() to split the input sequence (dna_seq) into triplets saved in s	      
        s = substr(dna_seq, i, 3)
       
        # make sure that input nt symbols are uppercase to match the hash keys
        s=toupper(s)

        # keep track of processed triplets (codons)
        triplet_counter++
       
        # check that the input corresponds to valid nucleotides
        if ( s ~ /[^acgtACGT]+/ ) { 
            print "WARNING: input triplet", triplet_counter, "=", s, 
        	  "contains a non-valid nucleotide symbol ..." > "/dev/stderr"
            break
        }

        # use the codon hash c as lookup table to translate the s triplet
        #   appending codons[s] to the growing peptide p
        { 
            # break out of loop if we get no more triplets 
            #	out of the input DNA string with substr()
            if (codons[s]=="") { 
               break
            }
            else if (s in codons == 0) { 
               # if the triplet is not contained in c, append "X" to p
               p = p unknown
            } 
            else { 
               # append aminoacid codons[s] to growing peptide
               p = p codons[s]
           }
        }
        i=i+3 # increment the counter of processed dna nucleotides by 3 
    }
    # run while loop until substring cannot retrieve any more triplets
    while (s!="")
        prots[header] = p
}
#---------------------------------------------------------------------------------------------------------

function print_help(prog, vers) 
{
  print "# AIMS:" > "/dev/stderr" 
  print "#  1. extracts the CDSs from a single GenBank file, saving them to genbank_basename_CDSs.fna" > "/dev/stderr"
  print "#  2. translates the CDSs and writes them to file genbank_basename_proteome.faa" > "/dev/stderr"
  print "#  3. writes the GenBank file in tabular format as genbank_basename.tsv" > "/dev/stderr"
  print "#  4. extracts the complete DNA string(s) from the input GenBank, saving them as genbank_basename.fsa" > "/dev/stderr"
  print "\n# NOTES:" > "/dev/stderr"
  print "#   1.", prog, "v"vers, "does not deal with CDSs containing introns" > "/dev/stderr"
  print "#        CDSs containing the join statement: complement(join(4497616..4498557,4498557..4498814)) are skipped" > "/dev/stderr"
  print "#        Use only for bacterial|mitochondrial|plastid|yeast genomes" > "/dev/stderr"
  print "\n# TODO:" > "/dev/stderr"
  print "#   1. add code to process split (pseudo)genes as in complement(join(4497616..4498557,4498557..4498814))" > "/dev/stderr"
  print "\n# USAGE of", prog, "v"vers":" > "/dev/stderr"
  printf "\t%s %s\n", prog, "genbank_file.gbk" > "/dev/stderr"  
  print "\tOR" > "/dev/stderr"  
  printf "\tfor f in *.gbk; do %s %s\n\n", prog, "$f; done" > "/dev/stderr"  
  
  exit 0
}
#---------------------------------------------------------------------------------------------------------
