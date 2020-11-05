#!/usr/bin/awk -f

# AUTHOR: Pablo Vinuesa, @pvinmex, https://www.ccg.unam.mx/~vinuesa/
# source: https://github.com/vinuesa/intro2linux
# translate_dna.awk VERSION:0.1
# AIM: translates a valid DNA string into proteins 
#      using the universal genetic code
# NOTE: this is a demo script to teach basic awk programming

BEGIN {
    
    progname = "translate_dna.awk"
    version  = 0.1  # nov 04, 2020
    
    if ( ARGC < 2 )
      Usage_Exit(progname, version)
    
    # initialize a hash named "c" holding the codon-aminoacid pairs 
    #   based on the universal genetic code
    c["ATA"]="I"; c["ATC"]="I"; c["ATT"]="I"; c["ATG"]="M";
    c["ACA"]="T"; c["ACC"]="T"; c["ACG"]="T"; c["ACT"]="T";
    c["AAC"]="N"; c["AAT"]="N"; c["AAA"]="K"; c["AAG"]="K";
    c["AGC"]="S"; c["AGT"]="S"; c["AGA"]="R"; c["AGG"]="R";
    c["CTA"]="L"; c["CTC"]="L"; c["CTG"]="L"; c["CTT"]="L";
    c["CCA"]="P"; c["CCC"]="P"; c["CCG"]="P"; c["CCT"]="P";
    c["CAC"]="H"; c["CAT"]="H"; c["CAA"]="Q"; c["CAG"]="Q";
    c["CGA"]="R"; c["CGC"]="R"; c["CGG"]="R"; c["CGT"]="R";
    c["GTA"]="V"; c["GTC"]="V"; c["GTG"]="V"; c["GTT"]="V";
    c["GCA"]="A"; c["GCC"]="A"; c["GCG"]="A"; c["GCT"]="A";
    c["GAC"]="D"; c["GAT"]="D"; c["GAA"]="E"; c["GAG"]="E";
    c["GGA"]="G"; c["GGC"]="G"; c["GGG"]="G"; c["GGT"]="G";
    c["TCA"]="S"; c["TCC"]="S"; c["TCG"]="S"; c["TCT"]="S";
    c["TTC"]="F"; c["TTT"]="F"; c["TTA"]="L"; c["TTG"]="L";
    c["TAC"]="Y"; c["TAT"]="Y"; c["TAA"]="*"; c["TAG"]="*";
    c["TGC"]="C"; c["TGT"]="C"; c["TGA"]="*"; c["TGG"]="W";
    
    unknown = "X   "  

 }
# -------------------- # 
# >>> MAIN PROGRAM <<< #  
# -------------------- # 
# Initialize variables: 
#  do-while loop control variable i (nt counter) 
#   and p, which will hold the translation product
{i=1; p=""; triplet_counter=0}

{
  # Here we run a do-while loop; the do loop is a variation of the while looping statement. 
  #  The do loop executes the body once and then repeats the body as long as the condition is true
  # We use the do-while loop, to get a first triplet string saved in s; 
  #  then the while loop keeps going until substr() got the last triplet, resulting in an empty s="".
  do {
         # First check that the script got some input
         #   if not, exit with an error message
         if(length($0) == 0) {
             print "ERROR: need a DNA sequence string to translate (valid DNA sequence, divisible by 3) "
             exit 1
         
         # Check that the DNA sequence string is divisible by 3 using the modulo operator
         #   if not, exit with an error message
         } else if(length($1)%3) { 
             print "ERROR: input DNA sequence not divisible by 3 ..."
             exit 2
         }

         # use substr() to split the input sequence ($1) into triplets saved in s         
         s=substr($1, i, 3)
         
         # keep track of processed triplets
         triplet_counter++
         
         # check that the input corresponds to valid nucleotides
         if ( s ~ /[^acgtACGT]+/ ) { 
             print "ERROR: input triplet", triplet_counter, "=", s, 
                      "contains a non-valid nucleotide symbol ..."
             exit 3
         }

         # make sure that input nt symbols are uppercase to match the hash keys
         s=toupper(s)
         
         # print the nucleotide sequence triplet, 
         #   followed by a space for easier visualization
         printf ("%s ", s)
         
         # use the codon hash c as lookup table to translate the s triplet
         #   appending c[s] to the growing peptide p
         { 
             # break out of loop if we get no more triplets 
             #   out of the input DNA string with substr()
             if (c[s]=="") { 
                break
             }
             else if (s in c == 0) { 
                # if the triplet is not contained in c, append "X   " to p
                p=p unknown
             } else { 
                # append aminoacid c[s] to growing peptide
                p=p c[s]"   "
             }
         }
         i=i+3 # increment the counter of processed dna nucleotides by 3 
    }
  # run while loop until substring cannot retrieve any more triplets
  while (s!="")
}

# this printf block prints the protein string
{ printf("\n %s\n", p) }


# function definition
function  Usage_Exit (prog, vers) # (prog, vers)
  {   
     print "USAGE:", prog, "v"vers
     print "   echo atggggtgttgtgggttgAAAGTGcccgggaaattaataCAG | ./translate_dna.awk -" > "/dev/stderr"
     print "   or: ./translate_dna.awk dna_string.txt" > "/dev/stderr"
     exit 1
  }


