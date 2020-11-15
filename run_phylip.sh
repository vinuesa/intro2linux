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

#: AIM: run PHYLIP's distance methods [NJ|UPGMA] for DNA and proteins with bootstrapping
#       This script was written to teach intermediate Bash scripting to my students at the 
#       Bachelor's Program in Genome Sciences at the Center for Genome Sciences, UNAM, Mexico
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
#	seqboot dnadist protdist neighbor consense nw_support
#   NOTE: Linux-compatible binaries are supplied in the distro\'s bin/ directory
#           https://github.com/vinuesa/intro2linux/tree/master/bin  

#: TODO:
#    implement also FM, ME, parsimony and ML analyses; 
#    check that the user provides proper PHYLIP files with at least 4 distinct sequences each

#: KNOWN BUGS: None.
#    Please report any errors you may encounter through the GitHub issue pages
#-------------------------------------------------------------------------------------------------------

# 0. Define strict bash settings
set -e              # exit on non-zero exit status
set -u              # exit if unset variables are encountered
set -o pipefail     # exit after unsuccessful UNIX pipe command

# set and export LC_NUMERIC=en_US.UTF-8, to avoid problems with locales tha use 1,32
LC_NUMERIC=en_US.UTF-8
export LC_NUMERIC

progname=${0##*/} # run_phylip.sh
VERSION=1.4 

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
    
    # v1.4 2020-11-14; * added function check_nw_utils_ok, to check that nw_display and nw_support can be executed
    #                    as they may be in path but cannot find /lib64/libc.so.6: version GLIBC_2.14
    #                    when the binaries are compiled with dynamic linking to the libs
    
    # v1.3 2020-11-14; * improved layout of output messages; 
    #                  * improved regex in extract_tree_from_outfile (now also for NJ tree)
    #                  * if nw_support and nw_display are not available, 
    #                     prints both the NJ and boot trees to screen if bootstrap analysis was performed
    
    # v1.2 2020-11-11; set and export LC_NUMERIC=en_US.UTF-8, to avoid problems with locales tha use 1,32
    
    # v1.1  2020-11-11; added function extract_tree_from_oufile and calls it on NJ|UPGMA bootRepl consensus trees
    #                      if nw_display is not available in PATH
    
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

function check_dependencies
{
    # check that the following PHYLIP binaries are in $PATH; die if not
    # optional dependencies found in bin/ dir: nw_support and nw_display
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
           grep GLIBC nw_utils_check.out.tmp &> /dev/null
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
    local matrix=$1
    
    # check that the distance matrix does not contain negative values due to wrong model parameters
    if grep -El ' \-[0-9\.]+' "$matrix"
    then
        echo "ERROR: computed negative distances!" >&2
	     [ "$runmode" -eq 1 ] && echo " You may need to ajust the model and|or gamma value; try lowering TiTv if > 6" >&2
	     [ "$runmode" -eq 2 ] && echo " You may need to ajust the matrix and|or gamma value" >&2
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
	     [ "$runmode" -eq 1 ] && echo "  You may need to adjust (lower?) TiTv=$TiTv" >&2
	     [ "$runmode" -eq 2 ] && echo "  You may need to use another matrix or adjunst gamma" >&2
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

function write_dnadist_params
{
    #get the input runmodes, models, bootstrap no. ...
    local model=$1
    local boot=$2
    local TiTv=$3
    local sequential=$4
    local gamma=$5
    local CV=$6

    # need to remove previous param file, since we are concatenating to it
    [ -s dnadist.params ] && rm dnadist.params
    
    #[ "$gamma" != 0 ] && CV=$( gamma2CV $gamma )

    # Runmode 1 = dnadist 
    [ "$model" = 'F84' ]          && touch dnadist.params
    [ -z "$model" ]               && touch dnadist.params 
    [ "$model" = 'Kimura' ]       && echo "D"                 >> dnadist.params
    [ "$model" = 'Jukes-Cantor' ] && echo -ne "D\nD\n"        >> dnadist.params
    [ "$model" = 'LogDet' ]       && echo -ne "D\nD\nD\n"     >> dnadist.params
    [ "$TiTv" != "2" ]            && echo -ne "T\n$TiTv\n"    >> dnadist.params
    [ "$gamma" != "0" ]           && echo -ne "G\n"           >> dnadist.params
    [ "$boot" -gt 0 ]             && echo -ne "M\nD\n$boot\n" >> dnadist.params
    [ "$sequential" -eq 1 ]       && echo -ne "I\n"           >> dnadist.params
                                     echo -ne "Y\n"           >> dnadist.params
    [ "$gamma" != "0" ]           && echo -ne "$CV\n"         >> dnadist.params
}
#---------------------------------------------------------------------------------#

function write_protdist_params
{
    #get the input runmodes, models and bootstrap no.
    # Models: JTT, PMB, PAM, Kimura
    local model=$1
    local boot=$2
    local sequential=$3
    local gamma=$4
    local CV=$5

    # need to remove previous param file, since we are concatenating to it
    [ -s protdist.params ] && rm protdist.params
    
    # Runmode 2 = protdist 
    [ "$model" = 'JTT' ]     && touch protdist.params
    [ -z "$model" ]          && touch protdist.params 
    [ "$model" = 'PMB' ]     && echo "P"                 >> protdist.params
    [ "$model" = 'PAM' ]     && echo -ne "P\nP\n"        >> protdist.params
    [ "$model" = 'Kimura' ]  && echo -ne "P\nP\nP\n"     >> protdist.params
    [ "$gamma" != "0" ]      && echo -ne "G\n"           >> protdist.params
    [ "$boot" -gt 0 ]        && echo -ne "M\nD\n$boot\n" >> protdist.params
    [ "$sequential" -eq 1 ]  && echo -ne "I\n"	         >> protdist.params
                                echo -ne "Y\n"           >> protdist.params
    [ "$gamma" != "0" ]      && echo -ne "$CV\n"         >> protdist.params				 

}
#---------------------------------------------------------------------------------#

function write_seqboot_params
{
    #get the input runmodes, models and bootstrap no.
    local boot=$1
    local sequential=$2
    
    # need to remove previous param file, since we are concatenating to it
    [ -s seqboot.params ] && rm seqboot.params

    # Write Seqboot params
    [ "$boot" -gt 0 ]       && echo -ne "R\n$boot\n" >> seqboot.params
    [ "$sequential" -eq 1 ] && echo -ne "I\n"        >> seqboot.params
                               echo -ne "Y\n"	     >> seqboot.params
			       echo -ne "$ROI\n"     >> seqboot.params   
}
#---------------------------------------------------------------------------------#

function write_neighbor_params
{
    #get the input runmodes, models and bootstrap no.
    local boot=$1
    local upgma=$2
    local outgroup=$3
    
    # need to remove previous param file, since we are concatenating to it
    [ -s neighbor.params ] && rm neighbor.params

    # Write $ROI params
    [ "$upgma" -gt 0 ]    && echo -ne "N\n"                 >> neighbor.params
    [ "$outgroup" -gt 1 ] && echo -ne "O\n$outgroup\n"      >> neighbor.params
    [ "$boot" -gt 0 ]     && echo -ne "M\n$boot\n$ROI\n"    >> neighbor.params
                             echo -ne "Y\n"                 >> neighbor.params
}
#---------------------------------------------------------------------------------#

function write_consense_params
{
    #very difficult
    [ -s consense.params ] && rm consense.params
    [ "$outgroup" -gt 1 ]  && echo -ne "O\n$outgroup\n"  >> consense.params 
                              echo -ne "Y\n"             >> consense.params
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
   #':b:c:m:I:o:R:t:hHDsuv'
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
   
 AIM: run PHYLIP's distance methods [NJ|UPGMA] for DNA and proteins with bootstrapping
      This code was written to teach basic Bash scripting to my students at the 
      Bachelor's Program in Genome Sciences, Center for Genome Sciences, UNAM, Mexico
      https://www.lcg.unam.mx/      

 OUTPUT: [NJ|UPGMA] phylogenies and, if requested, bootstrap consensus trees 
         and  [NJ|UPGMA] phylogenies with bootstrap support values mappend on bipartitions
 
 EXTERNAL DEPENDENCIES: 
   * PHYLIP (https://evolution.genetics.washington.edu/phylip.html) programs:
      seqboot dnadist protdist neighbor consense 
   * Newick utilities programs (http://cegg.unige.ch/newick_utils) programs:
      optional: nw_support and nw_display 
   - Note: Linux 64bit binaries available in the bin/ dir
  
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
while getopts ':b:g:i:m:R:o:t:hDHsuv' OPTIONS
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
   :)   printf "argument missing from -%s option\n" "$OPTARG"
   	 print_help
     	 exit 2 
     	 ;;
   \?)   echo "need the following args: "
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

# Check for correct model names


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
    gammaf=${gamma%.*}${gamma##*.}
else
    gammaf="$gamma"
fi    

# check model vs. runmode compatibility and suitable bootstrap values are provided
if [ "$runmode" -eq 1 ] && { [ "$model" = JTT ] || [ "$model" = PMG ] || [ "$model" = PAM ] || [ "$model" = "$def_prot_model" ]; }
then
       echo
       echo "# ERROR: $model is only valid when analyzing protein alignments under -R 2" >&2 
       echo
       print_help
       exit 1    
elif [ "$runmode" -eq 2 ] && { [ "$model" = F84 ] || [ "$model" = 'Jukes-Cantor' ] || [ "$model" = 'LogDet' ] || [ "$model" = "$def_DNA_model" ]; }
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

echo
echo "### $progname v.$VERSION run on $start_time with the following parameters:"
echo "# work_directory=$wkdir"
echo "# input_phylip=$input_phylip"
echo "# model=$model | gamma=$gamma | CV=$CV | gammaf=$gammaf | ti/tv ratio=$TiTv |
     outgroup=$outgroup | UPGMA=$upgma | bootstrap no.=$boot | ROI=$ROI"
echo "# runmode=$runmode"
echo "# DEBUG=$DEBUG" 
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


# 2.1) make sure we have an input file or die
if [ -s "$input_phylip" ]
then
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
    cp outfile "${input_phylip%.*}_${model}${gammaf}gamma_distMat.out" && \
      outfiles+=("${input_phylip%.*}_${model}${gammaf}gamma_distMat.out")
    
    # check that the distance matrix does not contain negative values due to wrong model parameters
    check_matrix "${input_phylip%.*}_${model}${gammaf}gamma_distMat.out"    
    
    mv outfile infile
elif [ "$runmode" -eq 2 ]
then
    echo "# running write_protdist_params $model 0 $sequential $gamma $CV"
    write_protdist_params "$model" "0" "$sequential" "$gamma" "$CV"
    echo "# running protdist < protdist.params"
    protdist < protdist.params &> /dev/null
    check_output outfile
    
    cp outfile "${input_phylip%.*}_${model}${gammaf}gamma_distMat.out" && \
      outfiles+=("${input_phylip%.*}_${model}${gammaf}gamma_distMat.out")
    
    # check that the distance matrix does not contain negative values due to wrong model parameters
    check_matrix "${input_phylip%.*}_${model}${gammaf}gamma_distMat.out"
    
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
     # nstead of appending one element, set the last item specifically, without any "unbound variable" error
     # t[${#t[*]}]=foo
     [ -s outtree ] && mv outtree "${input_phylip%.*}_${model}${gammaf}gamma_UPGMA.ph"	    && \
       outfiles[${#outfiles[*]}]=${input_phylip%.*}_${model}${gammaf}gamma_UPGMA.ph
     [ -s outfile ] && mv outfile "${input_phylip%.*}_${model}${gammaf}gamma_UPGMA.outfile" && \
       outfiles[${#outfiles[*]}]="${input_phylip%.*}_${model}${gammaf}gamma_UPGMA.outfile"
     
     upgma_tree="${input_phylip%.*}_${model}${gammaf}gamma_UPGMA.ph"
     echo "# wrote tree $upgma_tree to disk" 
       
     # check that there are no negative branch lengths in the nj_tree
     #  and display with nw_display, only if no bootstrapping is requested
     if [ "$boot" -eq 0 ] && [[ $(type -P nw_display) ]] && [[ "$nw_utils_ok" -eq 1 ]]
     then 
         display_treeOK "$upgma_tree"
     elif [ "$boot" -eq 0 ] && [[ ! $(type -P nw_display) ]]
     then
          echo "# extract_tree_from_outfile ${input_phylip%.*}_${model}${gammaf}gamma_UPGMA.outfile"
	  echo
	  extract_tree_from_outfile "${input_phylip%.*}_${model}${gammaf}gamma_UPGMA.outfile"
	  echo
     fi	 
else
     [ -s outtree ] && mv outtree "${input_phylip%.*}_${model}${gammaf}gamma_NJ.ph"	 && \
       outfiles[${#outfiles[*]}]="${input_phylip%.*}_${model}${gammaf}gamma_NJ.ph"
     [ -s outfile ] && mv outfile "${input_phylip%.*}_${model}${gammaf}gamma_NJ.outfile" && \
       outfiles[${#outfiles[*]}]="${input_phylip%.*}_${model}${gammaf}gamma_NJ.outfile"

     nj_tree="${input_phylip%.*}_${model}${gammaf}gamma_NJ.ph"
     echo "# wrote tree $nj_tree to disk"

     # check that there are no negative branch lengths in the nj_tree
     #  and display with nw_display, only if no bootstrapping is requested
     if [ "$boot" -eq 0 ] && [[ $(type -P nw_display) ]] && [[ "$nw_utils_ok" -eq 1 ]]
     then
         display_treeOK "$nj_tree"
     elif [ "$boot" -eq 0 ] && [[ ! $(type -P nw_display) ]]
     then
         echo "# extract_tree_from_outfile ${input_phylip%.*}_${model}${gammaf}gamma_NJ.outfile"
	 echo
	 extract_tree_from_outfile "${input_phylip%.*}_${model}${gammaf}gamma_NJ.outfile"
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
     
     if [ -s outtree ]
     then
	 # this is the file holding the trees for the n-distance matrices for n-boot replicated alignments
	 cp outtree "${input_phylip%.*}_${model}${gammaf}gamma_${boot}bootRepl_trees.nwk"
	 
	 # append to array with +=, otherwise will complain as unset with set -u 
	 # https://fvue.nl/wiki/Bash:_Error_%60Unbound_variable%27_when_appending_to_empty_array
	 outfiles+=("${input_phylip%.*}_${model}${gammaf}gamma_${boot}bootRepl_trees.nwk")
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
     
     if [ $upgma -gt 0 ]
     then
	 # https://fvue.nl/wiki/Bash:_Error_%60Unbound_variable%27_when_appending_to_empty_array
         # Set last item specifically
	 # nstead of appending one element, set the last item specifically, without any "unbound variable" error
	 # t[${#t[*]}]=foo
	 mv outtree "${input_phylip%.*}_UPGMAconsensus_${model}${gammaf}gamma_${boot}bootRepl.ph" && \
	   outfiles[${#outfiles[*]}]=${input_phylip%.*}_UPGMAconsensus_${model}${gammaf}gamma_${boot}bootRepl.ph
    	 mv outfile "${input_phylip%.*}_UPGMAconsensus_${model}${gammaf}gamma_${boot}bootRepl.outfile" && \
	  outfiles[${#outfiles[*]}]="${input_phylip%.*}_UPGMAconsensus_${model}${gammaf}gamma_${boot}bootRepl.outfile"
     else
	 mv outtree "${input_phylip%.*}_NJconsensus_${model}${gammaf}gamma_${boot}bootRepl.ph" && \
	  outfiles[${#outfiles[*]}]=${input_phylip%.*}_NJconsensus_${model}${gammaf}gamma_${boot}bootRepl.ph
    	 mv outfile "${input_phylip%.*}_NJconsensus_${model}${gammaf}gamma_${boot}bootRepl.outfile" && \
	  outfiles[${#outfiles[*]}]="${input_phylip%.*}_NJconsensus_${model}${gammaf}gamma_${boot}bootRepl.outfile"
     fi 


     # 5. Rename outtrees and tree outfiles; remap bootstrap values to bipartitions and display tree on screen
     if [ "$upgma" -gt 0 ]
     then
          # https://fvue.nl/wiki/Bash:_Error_%60Unbound_variable%27_when_appending_to_empty_array
          # Set last item specifically
          # nstead of appending one element, set the last item specifically, without any "unbound variable" error
          # t[${#t[*]}]=foo
          [ -s outtree ] && mv outtree "${input_phylip%.*}_${model}${gammaf}gamma_UPGMA.ph"	    && \
            outfiles[${#outfiles[*]}]=${input_phylip%.*}_${model}${gammaf}gamma_UPGMA.ph
          [ -s outfile ] && mv outfile "${input_phylip%.*}_${model}${gammaf}gamma_UPGMA.outfile" && \
            outfiles[${#outfiles[*]}]="${input_phylip%.*}_${model}${gammaf}gamma_UPGMA.outfile"
       
          # if we requested bootstrapping, map bootstrap values onto UPGMA tree using
          #   nw_support upgma.ph bootRepl_tree.ph > UPGMA_with_boot_support.ph
          if [ -s "${input_phylip%.*}_${model}${gammaf}gamma_${boot}bootRepl_trees.nwk" ] && [[ $(type -P nw_support) ]] && [ "$nw_utils_ok" -eq 1 ]
          then
	      echo "# mapping bootstrap values on UPGMA tree with nw_support ..."
	      nw_support "${input_phylip%.*}_${model}${gammaf}gamma_UPGMA.ph" \
	     "${input_phylip%.*}_${model}${gammaf}gamma_${boot}bootRepl_trees.nwk" > \
	     "${input_phylip%.*}_${model}${gammaf}gamma_UPGMA_with_${boot}boot_support.ph"

	      if [ -s "${input_phylip%.*}_${model}${gammaf}gamma_UPGMA_with_${boot}boot_support.ph" ] && [[ $(type -P nw_display) ]] && [ "$nw_utils_ok" -eq 1 ]
	      then
	          outfiles[${#outfiles[*]}]="${input_phylip%.*}_${model}${gammaf}gamma_UPGMA_with_${boot}boot_support.ph"
		  
		  nj_tree="${input_phylip%.*}_${model}${gammaf}gamma_UPGMA_with_${boot}boot_support.ph"
	          
		  # check that there are no negative branch lengths in the nj_tree 
		  #   before displaying with nw_display
		  display_treeOK "$nj_tree"
	      fi
	  elif [ -s "${input_phylip%.*}_UPGMAconsensus_${model}${gammaf}gamma_${boot}bootRepl.outfile" ] && { [[ ! $(type -P nw_support) ]] || [ "$nw_utils_ok" -eq 0 ]; }
	  then
               echo "# extract_tree_from_outfile ${input_phylip%.*}_${model}${gammaf}gamma_UPGMA.outfile"
	       echo
	       extract_tree_from_outfile "${input_phylip%.*}_${model}${gammaf}gamma_UPGMA.outfile"
	       echo
	       echo "# extract_tree_from_outfile ${input_phylip%.*}_UPGMAconsensus_${model}${gammaf}gamma_${boot}bootRepl.outfile"
	       echo
	       extract_tree_from_outfile "${input_phylip%.*}_UPGMAconsensus_${model}${gammaf}gamma_${boot}bootRepl.outfile"
	       echo
	  fi
     else
         [ -s outtree ] && mv outtree "${input_phylip%.*}_${model}${gammaf}gamma_NJ.ph"	 && \
           outfiles[${#outfiles[*]}]="${input_phylip%.*}_${model}${gammaf}gamma_NJ.ph"
         [ -s outfile ] && mv outfile "${input_phylip%.*}_${model}${gammaf}gamma_NJ.outfile" && \
           outfiles[${#outfiles[*]}]="${input_phylip%.*}_${model}${gammaf}gamma_NJ.outfile"

         # if we requested bootstrapping, map bootstrap values onto NJ tree using
         #   nw_support NJ.ph bootRepl_tree.ph > NJ_with_boot_support.ph
         if [ -s "${input_phylip%.*}_${model}${gammaf}gamma_${boot}bootRepl_trees.nwk" ] && [[ $(type -P nw_support) ]] && [ "$nw_utils_ok" -eq 1 ]
         then
             echo "# mapping bootstrap values on NJ tree with nw_support ..."
	     nw_support "${input_phylip%.*}_${model}${gammaf}gamma_NJ.ph" \
	     "${input_phylip%.*}_${model}${gammaf}gamma_${boot}bootRepl_trees.nwk" > \
	     "${input_phylip%.*}_${model}${gammaf}gamma_NJ_with_${boot}boot_support.ph"
	 
	     check_output "${input_phylip%.*}_${model}${gammaf}gamma_NJ_with_${boot}boot_support.ph"
	 
	     if [ -s "${input_phylip%.*}_${model}${gammaf}gamma_NJ_with_${boot}boot_support.ph" ] && [[ $(type -P nw_display) ]] && [ "$nw_utils_ok" -eq 1 ]
	     then
	         outfiles+=("${input_phylip%.*}_${model}${gammaf}gamma_NJ_with_${boot}boot_support.ph")

	         nj_tree="${input_phylip%.*}_${model}${gammaf}gamma_NJ_with_${boot}boot_support.ph"
	         
		 # check that there are no negative branch lengths in the nj_tree 
		 #   before displaying with nw_display
		 display_treeOK "$nj_tree"
	     fi
    	 elif [ -s "${input_phylip%.*}_NJconsensus_${model}${gammaf}gamma_${boot}bootRepl.outfile" ] && { [[ ! $(type -P nw_support) ]] || [ "$nw_utils_ok" -eq 0 ]; }
	 then

                echo "# extract_tree_from_outfile ${input_phylip%.*}_${model}${gammaf}gamma_NJ.outfile"
	        echo
	        extract_tree_from_outfile "${input_phylip%.*}_${model}${gammaf}gamma_NJ.outfile"
		echo

		echo "# extract_tree_from_outfile ${input_phylip%.*}_NJconsensus_${model}${gammaf}gamma_${boot}bootRepl.outfile"
		echo
		extract_tree_from_outfile "${input_phylip%.*}_NJconsensus_${model}${gammaf}gamma_${boot}bootRepl.outfile"
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
