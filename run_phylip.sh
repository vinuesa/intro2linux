#!/usr/bin/env bash

#-------------------------------------------------------------------------------------------------------
#: PROGRAM: run_phylip.sh
#: AUTHOR: Pablo Vinuesa, Center for Genomic Sciences, UNAM, Mexico
#:         https://www.ccg.unam.mx/~vinuesa/ twitter: @pvinmex
#
#: PROJECT START: October 16th, 2013
#    This program has been developed mainly for teaching purposes, 
#    with improvements/new features added as the script was used in 
#    diverse courses taught to undergrads at https://www.lcg.unam.mx
#    (LCG-UNAM) and the International Workshops on Bioinformatics (TIB)

#: AIM: run PHYLIP's distance methods [NJ|UPGMA] for DNA and proteins (dnadist|protdist) with optional bootstrapping
#       This script was written to teach intermediate Bash scripting to my students at the 
#       Bachelor's Program in Genome Sciences at the Center for Genome Sciences, UNAM, Mexico
#       https://www.lcg.unam.mx
#
#: INPUT: multiple sequence alignments (PROT|DNA) with at least 4 distinct sequences in phylip format
#
#: OUTPUT: [NJ|UPGMA] phylogenies and, if requested, bootstrap consensus trees

#: SOURCE: Freely available on GitHub @ https://github.com/vinuesa/intro2linux
#          Released under the GPLv3 License. 
#          http://www.gnu.org/copyleft/gpl.html

#### DEPENDENCIES 
# Assumes that the following binaries and scripts are all in $PATH, checked by check_dependencies()
#
#   1) Binaries from the PHYLIP package: 
#	seqboot dnadist protdist neighbor consense
#   NOTE: Linux-compatible PHYLIP binaries are supplied in the distro\'s bin/ directory
#           https://github.com/vinuesa/intro2linux/tree/master/bin  

#: TODO:
#    implement also parsimony and ML analyses and parallelize bootstrapping

#: KNOWN BUGS: None.
#    Please report any errors you may encounter through the GitHub issue pages
#-------------------------------------------------------------------------------------------------------

# make sure the user has at least bash version 4, since the script uses standard arrays (introduced in version 4),
#  but future development may require hashes, introduced in version 4
[ "${BASH_VERSION%%.*}" -lt 4 ] && echo "$HOSTNAME is running an ancient bash: ${BASH_VERSION}; ${0##*/} requires bash version 4 or higher" && exit 1

# 0. Define strict bash settings
set -e              # exit on non-zero exit status
set -u              # exit if unset variables are encountered
set -o pipefail     # exit after unsuccessful UNIX pipe command

# set and export LC_NUMERIC=en_US.UTF-8, to avoid problems with locales tha use 1,32
LC_NUMERIC=en_US.UTF-8
export LC_NUMERIC

args="$*"
progname=${0##*/} # run_phylip.sh
VERSION=2.0 

# GLOBALS
#DATEFORMAT_SHORT="%d%b%y" # 16Oct13
#TIMESTAMP_SHORT=$(date +${DATEFORMAT_SHORT})

date_F=$(date +%F |sed 's/-/_/g')-   # 2013_10_20
date_T=$(date +%T |sed 's/:/./g')    # 23.28.22 (hr.min.secs)
start_time="$date_F$date_T"
#current_year=$(date +%Y)

wkdir=$(pwd)

# initialize variables
def_DNA_model=F84
def_prot_model=JTT
input_phylip=
runmode=
boot=100
CV=1
model=
DEBUG=0
sequential=0
TiTv=2
upgma=0
outgroup=0
gamma="0.0"
outgroup=1

nw_utils_ok=

declare -a outfiles # to collect the output files written to disk

#---------------------------------------------------------------------------------#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION DEFINITIONS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
#---------------------------------------------------------------------------------#
function print_dev_history()
{
    cat <<EOF_
    
    $progname $VERSION has been developed mainly for teaching purposes, 
      with improvements/new features added as the script was used in 
      diverse courses taught to undergrads at https://www.lcg.unam.mx
      and the International Workshops on Bioinformatics (TIB)
      
    TODO: * implement also parsimony and ML analyses
          * parallelize bootstrapping  
    
    # v2.0 2020-11-19; * fixed bug in write_protdist_params in the if [ $model == JTT ] || [ ...] conditional
   		       * fixed an error in the model name input checks, changing PMG for PMB
   		       * fixed grep in check_nw_utils_ok and made more robust conditional checking to call the function
   		       * use variable expansion only once in naming outfiles outtress, saving it to a variable
   			 to make the code more streamlined and possibly a tiny bit faster
   		       * removed useless [ -s outtre|outfile ] tests for the orignal tress|outfiles in the bootstrapping block, 
   			   as these are already checked in the first code block computing tress from the original alignment
    
    # v1.6 2020-11-19; * streamlined write_PHYLIP_param functions by grouping echo calls in { } > params; avoiding concatenation 
    #                  * improved description of write_PHYLIP_param functionality

    # v1.5 2020-11-15; * added check_phylip_ok to validate input phylip file with -ge 4 sequences and -ge 4 characters
    		       * added remove_phylip_param_files
    		       * added option -I to call print_install_notes
    		       * added [ "${BASH_VERSION%%.*}" -lt 4 ] && die
    		       * makes more extensive use of internal Bash variables and variable expansions, including \${var/pattern/string}

    # v1.4 2020-11-14; * added function check_nw_utils_ok, to check that nw_display and nw_support can be executed
    			  as they may be in path but cannot find /lib64/libc.so.6: version GLIBC_2.14
    			  when the binaries are compiled with dynamic linking to the libs
    
    # v1.3 2020-11-14; * improved layout of output messages; 
    		       * improved regex in extract_tree_from_outfile (now also for NJ tree)
    		       * if nw_support and nw_display are not available, 
    			 prints both the NJ and boot trees to screen if bootstrap analysis was performed
    
    # v1.2 2020-11-11; set and export LC_NUMERIC=en_US.UTF-8, to avoid problems with locales that use 1,32 instead of 1.32
    
    # v1.1  2020-11-11; added function extract_tree_from_oufile and calls it on NJ|UPGMA bootRepl consensus trees
                          if nw_display is not available in PATH
    
    #v1.0  2020-11-09; This version has thorough error checking and reporting with improved code flow
         *  re-ingeneered the main code block. It now first runs the standard distance-matrix + clustering computations
	        before doing bootstrapping. This helps in rapidly detecting problematic model parameters, as reflected in negative distances
	 *  displays NJ|UPGMA tree without bootstrap values, if no bootstrapping is requested	
	 *  added functions:
	      - check_matrix to check if they contain negative values, issuing an error message
	      - display_treeOI display trees with nw_display, only if they do not contain negative branch lengths, issuing an error message
	      - print_dev_history to clean up the script\'s header section
	 *  made 100% shellcheck-compliant; prints all error messages to >&2 (STDERR)
         *  changed/fixed defaults for CV=1 & gamma="0.0"; 
	 *  use outfiles+=() to capture outfiles for -R 1 and -R 2 when boot -eq 0; now runs with -i -R1|2 as single opts
	 *  Changed default boot=100
	 *  now accepts TiTv=real like 6.7 and checks for real numbers passed to -t
	 *  checks if outtree contains negative branches before attempting to display with nw_display & prints proper message

    # v0.5 2020-11-09; added functions:
    	 *  rdm_odd_int to compute random odd integers to pass to PHYLIP programs
    	 *  check_output (silently checks for expected output files, dying if not found)

    #v0.4 2020-11-08; released to the public domain @ https://github.com/vinuesa/intro2linux
    	   * added nw_support and nw_display calls to map and display bottstrap support values onto NJ or UPGMA trees
    	   * fixed several warnings issued by shellcheck
    	   * added strict bash interpreter calling with set -e; set -u; set -o pipefail
    	   * localized variables received in functions
    	   * explicitly set default_DNA_model=F84 and default_prot_model=JTT
    	   * added option -v to print version and exit
    		
    #v0.3 September 15, 2017; improved checking of user input
    #v0.2 September 16, 2015; basic options and parameter checking for protdist phylogenies
    #v0.1 October 16, 2013; basic options and parameter checking for dnadist phylogenies
EOF_

exit 0
}
#---------------------------------------------------------------------------------#

function print_install_notes()
{
   cat <<INSTALL
   
    #1. PHYLIP PHYLogeny Inference Package by Joe Felensetein https://evolution.genetics.washington.edu/phylip.html
    
    The bin/ directory of the intro2linux distribution provides precompiled Linux-x86_64 binaries 
          for seqboot dnadist protdist neighbor consense. Copy them into your \$HOME/bin directory or any other in your PATH
    
    If you need executables for a different architecture, visit https://evolution.genetics.washington.edu/phylip/executables.html
    
    #2. Newick Utilities http://cegg.unige.ch/newick_utils
    The Newick Utilities have been described in an Open-Access paper (that is, available online for everyone):

    The Newick Utilities: High-throughput Phylogenetic tree Processing in the UNIX Shell
    Thomas Junier and Evgeny M. Zdobnov
    Bioinformatics 2010 26:1669-1670
    http://bioinformatics.oxfordjournals.org/content/26/13/1669.full
    doi:10.1093/bioinformatics/btq243
    
    The src/ dir contains the newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz source file
    
    Visit http://cegg.unige.ch/newick_utils if you need other pre-compiled versions
    
    Brief compilation notes
    
    1. Copy the newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz to a suitable directory in your \$HOME
    2. cd into the directory hoding the source *tar.gz file
    3. unpack and compile with the following commands 
          tar -xvzf newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz
	  cd newick-utils-1.6
	  ./configure --prefix=$HOME  # if you do NOT have administrator privileges
          #./configure                # if you have superuser privileges
	  make
	  make check
	  make install                # if you do NOT have administrator privileges
	  # sudo make install         # if you have superuser privileges

INSTALL

exit 0

}


function check_output()
{
    local outfile=$1
    if [ ! -s "$outfile" ]
    then
        echo
	echo " >>> ERROR! The expected output file $outfile was not produced, will exit now!" >&2
        echo
	exit 2
    fi
}
#---------------------------------------------------------------------------------#

function check_phylip_ok()
{
    # implements a rudimentary format validation check
    local phyfile=$1
    local num_tax=
    local num_char=

    num_tax=$(awk 'NR == 1 {print $1}' "$phyfile")
    num_char=$(awk 'NR == 1 {print $2}' "$phyfile")
    
    if [[ ! "$num_tax" =~ ^([0-9]+)$ ]] ||  [[ ! "$num_char" =~ ^([0-9]+)$ ]]
    then
         echo "ERROR: $phyfile is not a phylip-formatted input file!"
	 echo
	 exit 1  
    elif [[ "$num_tax" =~ ^([0-9]+)$ ]] && [ "$num_tax" -lt 4 ]
    then
         echo "ERROR: $phyfile should contain at least 4 taxa"
	 exit 1
    elif [[ "$num_char" =~ ^([0-9]+)$ ]] && [ "$num_char" -lt 5 ]
    then
         echo "ERROR: $phyfile should contain at least 5 aligned residues"
	 exit 1
    else
        echo "# File validation OK: $phyfile seems to be a standard phylip file ..."
        echo '-------------------------------------------------------------------------------------------'
    fi
}
#---------------------------------------------------------------------------------#

function check_dependencies
{
    # check that the following PHYLIP binaries are in $PATH; die if not
    # optional: the Newick_utilities src is found in the src/ dir
    dependencies=(seqboot dnadist protdist neighbor consense) 
    for programname in "${dependencies[@]}"
    do
       local bin=
       bin=$(type -P "$programname")
       if [ -z "$bin" ]; then
          echo
          echo "# ERROR: $programname not in place!" >&2
          echo "# ... you will need to install \"$programname\" first or include it in \$PATH" >&2
          echo "# ... exiting" >&2
          exit 1
       fi
    done

    echo
    echo '# Run check_dependencies() ... looks good: all required external dependencies are in place!'
    echo '-------------------------------------------------------------------------------------------'
}
#---------------------------------------------------------------------------------#

function check_nw_utils_ok()
{
    # This function checks that nw_display and nw_support can be executed
    #  as they may be in path but cannot find /lib64/libc.so.6: version GLIBC_2.14
    #  when the binaries are compiled with dynamic linking to the libs
    
    # The function returns 0 when nw_support nw_display cannot be executed and 1 when they can
    # The function does not run and return anything if nw_support nw_display are not in PATH
    
    dependencies=(nw_support nw_display)
    
    local nw_utils_ok=
    local bin=
    
    for programname in "${dependencies[@]}"
    do
        bin=$(type -P "$programname")
        if [ -n "$bin" ]; then
	   # prints something like nw_support: /lib64/libc.so.6: version `GLIBC_2.14' not found
	   "$programname" 2>> "nw_utils_check.out.tmp" 
        fi
    done
    
    # check the contents of nw_utils_check.out.tmp and set proper flag
    if [ -s nw_utils_check.out.tmp ]
    then
           # set the nw_utils_ok flag to 0 if the lib is not found, to check later in the code when nw_support and nw_display are called
           grep -i error nw_utils_check.out.tmp &> /dev/null
	   nw_utils_ok=$? 
	   [ "$nw_utils_ok" -eq 0 ] && echo "# WARNING: ${dependencies[*]} are in PATH but cannot find /lib64/libc.so.6: version GLIBC_2.14" >&2 
           rm nw_utils_check.out.tmp 
           echo  "$nw_utils_ok"
    fi
}
#---------------------------------------------------------------------------------#

function rdm_odd_int()
{
  # generate a random odd integer for seqboot and neighbor
  while true
  do 
     i=$RANDOM
     if (( "$i" % 2 )) # true when odd
     then 
         echo "$i" 
	 break
     fi
done

}
#---------------------------------------------------------------------------------#

function check_matrix()
{
    # check that the distance matrix does not contain negative values due to wrong model parameters
    local matrix=$1
    
    if grep -El ' \-[0-9\.]+' "$matrix"
    then
        echo "ERROR: computed negative distances!" >&2
	     [ "$runmode" -eq 1 ] && echo "You may need to ajust the model and|or gamma value; try lowering TiTv if > 6" >&2
	     [ "$runmode" -eq 2 ] && echo "You may need to ajust the matrix and|or gamma value" >&2
	exit 3
    fi
}
#---------------------------------------------------------------------------------#

function display_treeOK()
{
     # checks that tree does not contain branches with negative lengths, 
     #    as these cannot be displayed  with nw_dsiplay
     local tree=$1
     
     # check that there are no negative branch lengths in the nj_tree
     if ! grep -El '\-[0-9\.]+' "$tree"
     then 
     	 echo "# displaying $tree ..."
     	 echo
     	 nw_display -w 80 -Ir "$tree"
     	 echo
     else
     	 echo "ERROR: cannot display $nj_tree with nw_dsiplay, as it contains branches with negative lengths!" >&2
	     [ "$runmode" -eq 1 ] && echo "You may need to adjust (lower?) TiTv=$TiTv" >&2
	     [ "$runmode" -eq 2 ] && echo "You may need to use another matrix or adjunst gamma" >&2
     	 exit 4
     fi
}
#---------------------------------------------------------------------------------#

function extract_tree_from_outfile()
{
    # grep out lines containing tree characters | - and print to STDOUT
    local outfile=$1
   
    grep --color=never -E '[[:blank:]]+\||-' "$outfile" | grep -Ev 'Neighbor-|^---'
}

#------------------------------------- PHYLIP FUNCTIONS ---------------------------------------#
# >>> these are fucntions to write the command files to pass parameters to PHYLIP programs <<< # 

function remove_phylip_param_files()
{
   [ -s dnadist.params ] && rm dnadist.params
   [ -s protdist.params ] && rm protdist.params
   [ -s seqboot.params ] && rm seqboot.params
   [ -s neighbor.params ] && rm neighbor.params
   [ -s consense.params ] && rm consense.params
   return 0
}
#------------------------------------- PHYLIP FUNCTIONS ---------------------------------------#

function write_dnadist_params
{
    # writes a parameter file to run dnadist, based on provided arguments
    local model=$1
    local boot=$2
    local TiTv=$3
    local sequential=$4
    local gamma=$5
    local CV=$6

    # Runmode 1 = dnadist 
    if [ "$model" = 'F84' ] || [ -z "$model" ]
    then 
        {
          [ "$model" = 'Kimura' ]       && echo "D"                 
          [ "$model" = 'Jukes-Cantor' ] && echo -ne "D\nD\n"        
          [ "$model" = 'LogDet' ]       && echo -ne "D\nD\nD\n"     
          [ "$TiTv" != "2" ]            && echo -ne "T\n$TiTv\n"    
          [ "$gamma" != "0" ]           && echo -ne "G\n"           
          [ "$boot" -gt 0 ]             && echo -ne "M\nD\n$boot\n" 
          [ "$sequential" -eq 1 ]       && echo -ne "I\n"           
                                           echo -ne "Y\n"             
          [ "$gamma" != "0" ]           && echo -ne "$CV\n"         
    
        } > dnadist.params
   fi
}
#---------------------------------------------------------------------------------#

function write_protdist_params
{
    # write_protdist_params "$model" "$boot" "$sequential" "$gamma" "$CV"
    # writes a parameter file to run protdist, based on provided arguments
    # Models: JTT, PMB, PAM, Kimura
    local model=$1
    local boot=$2
    local sequential=$3
    local gamma=$4
    local CV=$5

    if [ "$model" = 'JTT' ] || [ -n "$model" ]
    then
       { # Runmode 2 = protdist 
         [ "$model" = 'PMB' ]     && echo "P"                 
         [ "$model" = 'PAM' ]     && echo -ne "P\nP\n"        
         [ "$model" = 'Kimura' ]  && echo -ne "P\nP\nP\n"     
         [ "$gamma" != "0" ]      && echo -ne "G\n"           
         [ "$boot" -gt 0 ]        && echo -ne "M\nD\n$boot\n" 
         [ "$sequential" -eq 1 ]  && echo -ne "I\n"    
                                     echo -ne "Y\n"       
         [ "$gamma" != "0" ]      && echo -ne "$CV\n"   			   
       } > protdist.params
    fi  
}
#---------------------------------------------------------------------------------#

function write_seqboot_params
{
    # writes a parameter file to run seqboot, based on provided arguments
    local boot=$1
    local sequential=$2
    
    # Write Seqboot params
    {
      [ "$boot" -gt 0 ]       && echo -ne "R\n$boot\n" 
      [ "$sequential" -eq 1 ] && echo -ne "I\n"        
                                 echo -ne "Y\n"
                                 echo -ne "$ROI\n"
    } > seqboot.params 			      
}
#---------------------------------------------------------------------------------#

function write_neighbor_params
{
    # writes a parameter file to run neighbor, based on provided arguments
    local boot=$1
    local upgma=$2
    local outgroup=$3
    
    {
       # Write $ROI params
       [ "$upgma" -gt 0 ]    && echo -ne "N\n"                 
       [ "$outgroup" -gt 1 ] && echo -ne "O\n$outgroup\n"      
       [ "$boot" -gt 0 ]     && echo -ne "M\n$boot\n$ROI\n"    
                                echo -ne "Y\n"        
    } > neighbor.params
}
#---------------------------------------------------------------------------------#

function write_consense_params
{
    #writes a parameter file to run consense; very difficult ;)
    [ -s consense.params ] && rm consense.params
    {
      [ "$outgroup" -gt 1 ]  && echo -ne "O\n$outgroup\n"  
                                echo -ne "Y\n"           
    } > consense.params				
}
#---------------------------------------------------------------------------------#

# don't leave litter behind ... remove intermediate input/output files
function cleanup_dir 
{
   [ -s infile ]  && rm infile
   [ -s outfile ] && rm outfile
   [ -s intree ]  && rm intree
   [ -s outtree ] && rm outtree
   
   for file in *.params
   do
       [ -s "$file" ] && rm "$file"
   done    
}

#-------------------------------- END PHYLIP FUNCTIONS----------------------------#

function print_help
{
   #':b:c:m:I:o:R:t:hHIDsuv'
   cat<<EOF
   $progname v.$VERSION [OPTIONS]
   
   -h  prints this help message
   
 REQUIRED
   -i <input_phylip_file> (if sequential, use also flag -s)
   -R <integer>  RUNMODE
        1    run dnadist  (expects DNA  alignment in phylip format)
	2    run protdist (expects PROT alignment in phylip format)

 OPTIONAL
   -b <integer> No. of bootstrap pseudoreplicates [default: $boot]
   -g <real number> alpha for gamma distribution  [default: $gamma]
   -o <integer> outgroup sequence no.             [default: $outgroup]
   -s sequential format; -s (flag; no val)!       [default: interleaved]
   -t <digit> Transition/Transversion ratio       [default: $TiTv]
   -m <model name>  NOTE: TYPE AS SHOWN!
        DNA:  F84, Kimura, Jukes-Cantor, LogDet   [default: $def_DNA_model]
        PROT: JTT, PMB, PAM, Kimura               [default: $def_prot_model]
   -u <flag> use UPGMA for clustering             [default:  NJ]
   -v <flag> print program version and exit
   -D <flag> Activate debugging to keep cmd files [default: $DEBUG]
   -H <flag> print development history and exit   
   -I <flag> print installation notes and exit   
   
 AIM: run PHYLIP\'s distance methods [NJ|UPGMA] for DNA and proteins with bootstrapping
      This code was written to teach basic Bash scripting to my students at the 
      Bachelor\'s Program in Genome Sciences, Center for Genome Sciences, UNAM, Mexico
      https://www.lcg.unam.mx/      

 OUTPUT: [NJ|UPGMA] phylogenies and, if requested, bootstrap consensus trees 
         and  [NJ|UPGMA] phylogenies with bootstrap support values mappend on bipartitions
 
 EXTERNAL DEPENDENCIES: 
   * PHYLIP (https://evolution.genetics.washington.edu/phylip.html) programs:
      seqboot dnadist protdist neighbor consense 
   * Newick utilities programs (http://cegg.unige.ch/newick_utils) programs:
      optional: nw_support and nw_display; need to install separately
   - Notes: PHYLIP Linux 64bit binaries available in the bin/ dir
            Copy them to your \$HOME/bin directory or any other in your PATH
  
 LICENSING & SOURCE
       Author: Pablo Vinuesa | https://www.ccg.unam.mx/~vinuesa/ | twitter: @pvinmex
       Released under the GNU General Public License version 3 (GPLv3)
       http://www.gnu.org/copyleft/gpl.html
       source: https://github.com/vinuesa/intro2linux

EOF

exit 1   
}

#-------------------------------------------------------------------------------------------------#
#---------------------------------------- GET OPTIONS --------------------------------------------#
#-------------------------------------------------------------------------------------------------#

# GETOPTS
while getopts ':b:g:i:m:R:o:t:hDHIsuv' OPTIONS
do
   case $OPTIONS in

   b)   boot=$OPTARG
        ;;
   g)   gamma=$OPTARG
	[ "$gamma" != 0 ] && CV=$(echo "1/sqrt($gamma)" | bc -l) && printf -v CV "%.4f\n" "$CV"
        ;;
   h)   print_help
        ;;
   i)   input_phylip=$OPTARG
        ;;
   s)   sequential=1
        ;;
   m)   model=$OPTARG
        ;;
   o)   outgroup=$OPTARG
        ;;
   R)   runmode=$OPTARG
        ;;
   t)   TiTv=$OPTARG
        ;;
   u)   upgma=1
        ;;
   v)   echo "$progname v$VERSION" && exit
        ;;
   D)   DEBUG=1
        ;;
   H)   print_dev_history
        ;;
   I)   print_install_notes
        ;;
   :)   printf "argument missing from -%s option\n" "$OPTARG"
   	print_help
     	exit 2 
     	;;
   \?)  echo "need the following args: "
   	print_help
        exit 3
	;;
   esac >&2   # print the ERROR MESSAGES to STDERR

done

shift $((OPTIND - 1))

#--------------------------#
# >>> Check User Input <<< #
#--------------------------#

if [ -z "$input_phylip" ]
then
       echo
       echo "# ERROR: no input phylip file defined!"
       print_help
       exit 1    
fi

if [ -z "$runmode" ]
then
       echo
       echo "# ERROR: no runmode defined!"
       print_help
       exit 1    
fi

# automatically set TiTv=0 when running with protein matrices or the LogDet DNA model
[ "$runmode" -eq 2 ] && TiTv=0
[ "$model" == LogDet ] && TiTv=0


# check that bootstrap value is an integer
re='^[0-9]+$'
if [[ ! "$boot" =~ $re ]]
then
   echo
   echo "# ERROR: boot:$boot is not a positive integer >= 0; provide a value between 100 and 1000" >&2 
   echo
   print_help
   echo
   exit 3
fi

# check that Ti/Tv is a real number, integer or decimal
re_TiTv='^[0-9.]+$'
if [[ ! "$TiTv" =~ $re_TiTv ]] && [ "$runmode" -eq 1 ]
then
   echo
   echo "# ERROR: TiTv:$TiTv is not an integer >= 0; provide a value between 0-10" >&2 
   echo
   print_help
   echo
   exit 3
elif [[ "$TiTv" =~ $re_TiTv ]] && [ "$runmode" -eq 1 ] && { [ "$model" == Jukes-Cantor ] || [ "$model" == LogDet ]; }
then
       echo
       echo "# ERROR: $model is only valid when analyzing DNA alignments under the F84 or K2P models" >&2
       echo
       print_help
       exit 1    
fi

if [ "$gamma" != 0 ] # note, need to treat as string, since gamma will be generalle a float like 0.5
then
    # this is to avoid having dots within filename (converts gamma==0.2 to gammaf=02)
    gammaf=${gamma/\./} #gammaf=${gamma%.*}${gamma##*.}
else
    gammaf="$gamma"
fi    

# check model vs. runmode compatibility and suitable bootstrap values are provided
#   using two alternative sintaxes for educational purposes
if [ "$runmode" -eq 1 ] && { [ "$model" = JTT ] || [ "$model" = PMB ] || [ "$model" = PAM ] || [ "$model" = "$def_prot_model" ]; }
then
       echo
       echo "# ERROR: $model is only valid when analyzing protein alignments under -R 2" >&2 
       echo
       print_help
       exit 1    
elif [ "$runmode" -eq 2 ] && [[ "$model" =~ ^(F84|Jukes-Cantor|LogDet|"$def_DNA_model")$ ]]
then
       echo
       echo "# ERROR: $model is only valid when analyzing DNA alignments under -R 1" >&2
       echo
       print_help
       exit 1    
elif [ "$boot" -lt 0 ]
then
       echo
       echo "# ERROR: bootstrap value=$boot is < 0 and not permitted" >&2
       echo
       print_help
       exit 1    
elif [ "$boot" -gt 1000 ]
then
       echo
       echo "# WARNING: bootstrap value=$boot is > 1000. This may take a long time to run."
       echo -n "# Are you sure you want to continue anyway? Type Y|N "
       read -r answer
       if [ "$answer" = N ] || [ "$answer" = n ] 
       then
           echo
	   echo "# Ok, will exit then!"
	   echo
	   exit 2
       else
           echo
	   echo "# Ok, will proceed running with $boot bootstrap pseudoreplicates ..."
	   echo
       fi 
elif [ "$boot" -gt 0 ] && [ "$boot" -lt 100 ]
then
       echo
       echo "# WARNING: bootstrap value=$boot is a bit low. Use a value >= 100."
       echo -n "# Are you sure you want to continue anyway? Type Y|N "
       read -r answer
       if [ "$answer" = N ] || [ "$answer" = n ] 
       then
           echo
	   echo "# Ok, will exit then!"
	   echo
	   exit 2
       else
           echo
	   echo "# Ok, will proceed ..."
	   echo
       fi 
fi

# set the default DNA or PROT models, if not defined by the user
if [ "$runmode" -eq 1 ] && [ -z "$model" ]
then
       model="$def_DNA_model"
       #echo "# DNA model set to default: $def_DNA_model..."       
fi

if [ "$runmode" -eq 2 ] && [ -z "$model" ]

then
       model="$def_prot_model"
       #echo "# Protein matrix set to default: $def_prot_model ..."       
fi

# check that the user provided a valid DNA substitution model
if [ "$runmode" -eq 1 ] && [[ ! "$model" =~ ^(F84|Kimura|Jukes-Cantor|LogDet|"$def_DNA_model")$ ]]
then
   echo
   echo "# ERROR: $model is not a recognized substitution model for DNA sequences used by PHYLIP" >&2 
   echo
   print_help
   echo
   exit 3
fi

# check that the user provided a valid name for empirical substitution matrix
#    note the much shorter and cleaner notation of this test using extended regexes, than the previous one
if [ "$runmode" -eq 2 ] && [[ ! "$model" =~ ^(JTT|PMB|PAM|Kimura)$ ]]
then
   echo
   echo "# ERROR: $model is not a recognized substitution matrix for protein sequences used by PHYLIP" >&2 
   echo
   print_help
   echo
   exit 3
fi

#>>> Set the script's run random odd number required by diverse PHYLIP programs <<<
ROI=$(rdm_odd_int)

# print the run settings
echo
echo "### $progname v.$VERSION run on $start_time with the following parameters:"
echo "# work_directory=$wkdir"
echo "# input_phylip=$input_phylip"
echo "# model=$model | gamma=$gamma | CV=$CV | gammaf=$gammaf | ti/tv ratio=$TiTv |
     outgroup=$outgroup | UPGMA=$upgma | bootstrap no.=$boot | ROI=$ROI"
echo "# runmode=$runmode"
echo "# DEBUG=$DEBUG" 
echo "# command: $progname ${args[*]}"
echo

#-------------------------------------------------------------------------------------------------#
#------------------------------------------ MAIN CODE --------------------------------------------#
#-------------------------------------------------------------------------------------------------#

# 1. make sure the external dependencies are found in PATH
#     and that nw_display and nw_support find the required GLIBC_2.14 in /lib64/libc.so.6
check_dependencies
nw_utils_ok=$(check_nw_utils_ok) # nw_utils_ok -eq 1 when OK, -eq 0 when not

# 2. Start processing the input file
# make sure there are no old outfile or outtree files lying around from previous runs
[ -s infile ] && rm infile
[ -s outfile ] && rm outfile
[ -s outtree ] && rm outtree
# nor old params files
remove_phylip_param_files

# 2.1) make sure we have an input file or die
if [ -s "$input_phylip" ]
then
     # make basic format validation check
     check_phylip_ok "$input_phylip"
     cp "$input_phylip" infile
else
     echo
     echo "# FATAL ERROR: input phylip file $input_phylip does not exist or is empty" >&2 
     echo "# Make sure $input_phylip is in $wkdir. Exiting ..." >&2
     echo
     exit 1
fi

# ----------------------------------------------------------- #
# >>>>>>>>>>>>>>>> Compute distance matrices <<<<<<<<<<<<<<<< #
# ----------------------------------------------------------- #

# 3. In any case (with or without bootstrap) compute the NJ tree for the original phylip file
#     Run dnadist or protdist, as required   
echo
echo ">>> Computing distance matrix for $input_phylip ..."

if [ "$runmode" -eq 1 ]		 
then
    echo "# running write_dnadist_params $model 0 $TiTv $sequential $gamma $CV"
    write_dnadist_params "$model" "0" "$TiTv" "$sequential" "$gamma" "$CV"
    echo "# running dnadist < dnadist.params"
    dnadist < dnadist.params &> /dev/null
    check_output outfile
    
    # https://fvue.nl/wiki/Bash:_Error_%60Unbound_variable%27_when_appending_to_empty_array
    # Set last item specifically
    # nstead of appending one element, set the last item specifically, without any "unbound variable" error
    # t[${#t[*]}]=foo
    dnadist_outfile="${input_phylip%.*}_${model}${gammaf}gamma_distMat.out"
    cp outfile "$dnadist_outfile" && \
      outfiles+=("$dnadist_outfile")
    
    # check that the distance matrix does not contain negative values due to wrong model parameters
    check_matrix "$dnadist_outfile"    
    
    mv outfile infile
elif [ "$runmode" -eq 2 ]
then
    echo "# running write_protdist_params $model 0 $sequential $gamma $CV"
    write_protdist_params "$model" "0" "$sequential" "$gamma" "$CV"
    echo "# running protdist < protdist.params"
    protdist < protdist.params &> /dev/null
    check_output outfile
    
    protdist_outfile="${input_phylip%.*}_${model}${gammaf}gamma_distMat.out"
    cp outfile "$protdist_outfile" && \
      outfiles+=("$protdist_outfile")
    
    # check that the distance matrix does not contain negative values due to wrong model parameters
    check_matrix "$protdist_outfile"
    
    mv outfile infile
fi


# ------------------------------------------------------------------------------------ #
# >>>>>>>>>>>>>>>> Computing NJ|UPGMA trees from original alignments <<<<<<<<<<<<<<<<< #
# ------------------------------------------------------------------------------------ #

echo
echo ">>> Computing distance tree for $input_phylip ..."  

# 4. now that we have the dist matrix, do the clustering with NJ or UPGMA
echo "# running write_neighbor_params 0 $upgma"
write_neighbor_params "0" "$upgma" "$outgroup"
echo "# running neighbor < neighbor.params"
neighbor < neighbor.params &> /dev/null
check_output outfile
check_output outtree

# 5.1 rename outtrees and tree outfiles; remap bootstrap values to bipartitions and display tree to screen
if [ "$upgma" -gt 0 ]
then
     # https://fvue.nl/wiki/Bash:_Error_%60Unbound_variable%27_when_appending_to_empty_array
     # Set last item specifically
     # instead of appending one element, set the last item specifically, without any "unbound variable" error
     # t[${#t[*]}]=foo
     upgma_tree=
     upgma_outfile=
     
     if [ -s outtree ] 
     then
          upgma_tree="${input_phylip%.*}_${model}${gammaf}gamma_UPGMA.ph"
          mv outtree "$upgma_tree"
	  echo "# wrote tree $upgma_tree to disk" 
	  outfiles[${#outfiles[*]}]="$upgma_tree"
     fi 
     
     if [ -s outfile ]
     then
         upgma_outfile="${input_phylip%.*}_${model}${gammaf}gamma_UPGMA.outfile"
         mv outfile "$upgma_outfile"
         outfiles[${#outfiles[*]}]="$upgma_outfile"
     fi 
       
     # check that there are no negative branch lengths in the nj_tree
     #  and display with nw_display, only if no bootstrapping is requested
     if [ "$boot" -eq 0 ] && [[ $(type -P nw_display) ]] && [[ "$nw_utils_ok" -eq 1 ]]
     then 
         display_treeOK "$upgma_tree"
     elif [ "$boot" -eq 0 ] && { [[ ! $(type -P nw_display) ]] || [[ "$nw_utils_ok" -ne 1 ]]; }
     then
          echo "# extract_tree_from_outfile $upgma_outfile"
	  echo
	  extract_tree_from_outfile "$upgma_outfile"
	  echo
     fi	 
else
     nj_tree=
     nj_outfile=

     if [ -s outtree ] 
     then
	  nj_tree="${input_phylip%.*}_${model}${gammaf}gamma_NJ.ph"
          mv outtree "$nj_tree"
          echo "# wrote tree $nj_tree to disk"
          outfiles[${#outfiles[*]}]="$nj_tree"
     fi
     
     if [ -s outfile ] 
     then
          nj_outfile="${input_phylip%.*}_${model}${gammaf}gamma_NJ.outfile"
          mv outfile "$nj_outfile"
          outfiles[${#outfiles[*]}]="$nj_outfile"
     fi

     # check that there are no negative branch lengths in the nj_tree
     #  and display with nw_display, only if no bootstrapping is requested
     if [ "$boot" -eq 0 ] && [[ $(type -P nw_display) ]] && [[ "$nw_utils_ok" -eq 1 ]]
     then
         display_treeOK "$nj_tree"
     elif [ "$boot" -eq 0 ] && { [[ ! $(type -P nw_display) ]] || [[ "$nw_utils_ok" -ne 1 ]]; }
     then
         echo "# extract_tree_from_outfile $nj_outfile"
	 echo
	 extract_tree_from_outfile "$nj_outfile"
	 echo
     fi	 
fi

echo "# > finished computing distance matrix and tree for $input_phylip!"

if [ "$boot" -eq 0 ]
then
    # 6. Print final output summary message
    echo
    echo '===================== OUTPUT SUMMARY ====================='
    no_outfiles=${#outfiles[@]}
    echo "# $no_outfiles output files were generated:"
    printf "%s\n" "${outfiles[@]}"

    echo
    echo -n "# FINISHED run at: "; date 
    echo "   ==> Exiting now ..."
    echo
else
    echo "=================================================================="
    echo
fi


# ----------------------------------------------------------- #
# >>>>>>>>>>>>>>>>>>>> Bootstrap Analysis <<<<<<<<<<<<<<<<<<< #
# ----------------------------------------------------------- #

# Run seqboot if requested
# run seqboot if -b > 0

if [ "$boot" -gt 0 ] 
then
     # 1. restore the original infile for the bootstrapping and standard NJ/UPGMA analysis below
     cp "$input_phylip" infile
     check_output infile

     echo
     echo ">>> Bootstrap Analysis based on $boot pseudoreplicates for $input_phylip ..."
     echo "# running  write_seqboot_params $boot $sequential"
     write_seqboot_params "$boot" "$sequential"
     echo "# running seqboot < seqboot.params &> /dev/null"
     seqboot < seqboot.params &> /dev/null
     check_output outfile
     mv outfile infile

     echo "# > computing distance matrices on $boot bootstrapped alignments ..."
     # 2. if bootstrapping, then compute consensus tree
     # we need to run dnadist or protdist, depending on runmode
     if [ "$runmode" -eq 1 ]
     then
	 echo "# running write_dnadist_params $model $boot $TiTv $sequential $gamma $CV"
    	 write_dnadist_params "$model" "$boot" "$TiTv" "$sequential" "$gamma" "$CV"
	 echo "# running dnadist < dnadist.params &> /dev/null"
    	 dnadist < dnadist.params &> /dev/null
	 check_output outfile
	 
	 # check that matrix does not contain negative values
	 check_matrix outfile
	 
	 mv outfile infile
     elif [ "$runmode" -eq 2 ]
     then
	 echo
	 echo "# running write_protdist_params $model $boot $sequential $gamma $CV"
    	 write_protdist_params "$model" "$boot" "$sequential" "$gamma" "$CV"
	 echo "# running protdist < protdist.params &> /dev/null"
    	 protdist < protdist.params &> /dev/null
	 check_output outfile
	 
	 # check that matrix does not contain negative values
	 check_matrix outfile

	 mv outfile infile
     fi

     # 3. Now we have the distance matrices and we can proceed equaly for runmodes 1 and 2
     #    >>> Run neighbor
     echo "# > Computing distance trees from bootstrapped data ..."	  
     echo "# running write_neighbor_params $boot $upgma $outgroup"
     write_neighbor_params "$boot" "$upgma" "$outgroup"
     echo "# running neighbor < neighbor.params &> /dev/null"
     neighbor < neighbor.params &> /dev/null
     check_output outtree
     
     boot_trees=
     if [ -s outtree ]
     then
	 # this is the file holding the trees for the n-distance matrices for n-boot replicated alignments
	 boot_trees="${input_phylip%.*}_${model}${gammaf}gamma_${boot}bootRepl_trees.nwk"
	 cp outtree "$boot_trees"
	 
	 # append to array with +=, otherwise will complain as unset with set -u 
	 # https://fvue.nl/wiki/Bash:_Error_%60Unbound_variable%27_when_appending_to_empty_array
	 outfiles+=("$boot_trees")
     fi
     mv outtree intree
     rm outfile

     # 4. Compute consensus tree with consense
     echo "# > Computing MJR consensus tree from trees reconstructed from bootstrap pseudoreplicates ..."
     echo "# running write_consense_params"
     write_consense_params
     echo "# running consense < consense.params &> /dev/null"
     consense < consense.params &> /dev/null
     check_output outtree
     check_output outfile
     
     # variables holding consensus trees
     upgma_consensus_tree=
     upgma_consensus_outfile=
     nj_consensus_tree=
     nj_consensus_outfile=

     if [ "$upgma" -gt 0 ]
     then
	 # https://fvue.nl/wiki/Bash:_Error_%60Unbound_variable%27_when_appending_to_empty_array
         # Set last item specifically
	 # instead of appending one element, set the last item specifically, without any "unbound variable" error
	 # t[${#t[*]}]=foo
	 upgma_consensus_tree="${input_phylip%.*}_UPGMAconsensus_${model}${gammaf}gamma_${boot}bootRepl.ph"
	 mv outtree "$upgma_consensus_tree"
	 outfiles[${#outfiles[*]}]="$upgma_consensus_tree"
    	 
	 upgma_consensus_outfile="${input_phylip%.*}_UPGMAconsensus_${model}${gammaf}gamma_${boot}bootRepl.outfile"
	 mv outfile "$upgma_consensus_outfile"
	 outfiles[${#outfiles[*]}]="$upgma_consensus_outfile"
     else
         nj_consensus_tree="${input_phylip%.*}_NJconsensus_${model}${gammaf}gamma_${boot}bootRepl.ph"
	 mv outtree "$nj_consensus_tree"
	 outfiles[${#outfiles[*]}]="$nj_consensus_tree"
    	 
	 nj_consensus_outfile="${input_phylip%.*}_NJconsensus_${model}${gammaf}gamma_${boot}bootRepl.outfile"
	 mv outfile "$nj_consensus_outfile"
	 outfiles[${#outfiles[*]}]="$nj_consensus_outfile"
     fi 

     # 5. Rename outtrees and tree outfiles; remap bootstrap values to bipartitions and display tree on screen
     if [ "$upgma" -gt 0 ]
     then
          # if we requested bootstrapping, map bootstrap values onto UPGMA tree using
          #   nw_support upgma.ph bootRepl_tree.ph > UPGMA_with_boot_support.ph
          if [ -s "$upgma_tree" ] && [ -s "$boot_trees" ] && [[ $(type -P nw_support) ]] && [ "$nw_utils_ok" -eq 1 ]
          then
	      upgma_tree_with_boot="${input_phylip%.*}_${model}${gammaf}gamma_UPGMA_with_${boot}boot_support.ph"
	      echo "# mapping bootstrap values on UPGMA tree with nw_support ..."
	      nw_support "$upgma_tree" "$boot_trees" > "$upgma_tree_with_boot"

	      if [ -s "$upgma_tree_with_boot" ] && [[ $(type -P nw_display) ]] && [ "$nw_utils_ok" -eq 1 ]
	      then
	          outfiles[${#outfiles[*]}]="$upgma_tree_with_boot"
	          
		  # check that there are no negative branch lengths in the nj_tree 
		  #   before displaying with nw_display
		  display_treeOK "$upgma_tree_with_boot"
	      fi
	  elif [ -s "$upgma_outfile" ] && [ -s "$upgma_consensus_outfile" ] && { [[ ! $(type -P nw_support) ]] || [ "$nw_utils_ok" -ne 1 ]; }
	  then
               echo "# extract_tree_from_outfile $upgma_outfile"
	       echo
	       extract_tree_from_outfile "$upgma_outfile"
	       echo
	       echo "# extract_tree_from_outfile $upgma_consensus_outfile"
	       echo
	       extract_tree_from_outfile "$upgma_consensus_outfile"
	       echo
	  fi
     else
         # if we requested bootstrapping, map bootstrap values onto NJ tree using
         #   nw_support NJ.ph bootRepl_tree.ph > NJ_with_boot_support.ph
         if [ -s "$nj_tree" ] && [ -s "$boot_trees" ] && [[ $(type -P nw_support) ]] && [ "$nw_utils_ok" -eq 1 ]
         then
             nj_tree_with_boot="${input_phylip%.*}_${model}${gammaf}gamma_NJ_with_${boot}boot_support.ph"
	     echo "# mapping bootstrap values on NJ tree with nw_support ..."
	     nw_support "$nj_tree" "$boot_trees" > "$nj_tree_with_boot"
	 
	     check_output "$nj_tree_with_boot"
	 
	     if [ -s "$nj_tree_with_boot" ] && [[ $(type -P nw_display) ]] && [ "$nw_utils_ok" -eq 1 ]
	     then
	         outfiles+=("$nj_tree_with_boot")
	         
		 # check that there are no negative branch lengths in the nj_tree 
		 #   before displaying with nw_display
		 display_treeOK "$nj_tree_with_boot"
	     fi
    	 elif [ -s "$nj_outfile" ] && [ -s "$nj_consensus_outfile" ] && { [[ ! $(type -P nw_support) ]] || [ "$nw_utils_ok" -ne 1 ]; }
	 then

                echo "# extract_tree_from_outfile $nj_outfile"
	        echo
	        extract_tree_from_outfile "$nj_outfile"
		echo

		echo "# extract_tree_from_outfile $nj_consensus_outfile"
		echo
		extract_tree_from_outfile "$nj_consensus_outfile"
		echo
         fi
    fi
fi

# 6. Tidy up: remove the *params files and other temporary files 
# that could interfere with future runs and litter the directory

[ "$DEBUG" -eq 0 ] && cleanup_dir

# 7. Print final output summary message
echo
echo '===================== OUTPUT SUMMARY ====================='

no_outfiles=${#outfiles[@]}
echo "# $no_outfiles output files were generated:"
printf "%s\n" "${outfiles[@]}"

echo
echo -n "# FINISHED run at: "; date 
echo "   ==> Exiting now ..."
echo
