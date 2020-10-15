#!/usr/bin/env bash

#: PROGRAM: XXX
#: AUTHOR: Pablo Vinuesa, Center for Genomic Sciences, UNAM, Mexico
#:         http://www.ccg.unam.mx/~vinuesa/
#: LICENSE: GNU General Public License v3.0
#  Permissions of this strong copyleft license are conditioned on making available complete source code of licensed works and modifications, which include larger works using a licensed work, under the same license. Copyright and license notices must be preserved. Contributors provide an express grant of patent rights.

# GLOBALS
args="$@"
progname=${0##*/} # XXX
VERSION=0.1 # XXX

# TODO

wkdir=$(pwd)

DATEFORMAT_SHORT="%d%b%y"
TIMESTAMP_SHORT=$(date +${DATEFORMAT_SHORT})

date_F=$(date +%F |sed 's/-/_/g')-
date_T=$(date +%T |sed 's/:/./g')
start_time="$date_F$date_T"
current_year=$(date +%Y)

#>>> set color in bash 
#  SEE: echo http://stackoverflow.com/questions/5947742/how-to-change-the-output-color-of-echo-in-linux
# ANSI escape codes
# Black        0;30     Dark Gray     1;30
# Red          0;31     Light Red     1;31
# Green        0;32     Light Green   1;32
# Brown/Orange 0;33     Yellow        1;33
# Blue         0;34     Light Blue    1;34
# Purple       0;35     Light Purple  1;35
# Cyan         0;36     Light Cyan    1;36
# Light Gray   0;37     White         1;37

RED='\033[0;31m'
GREEN='\033[0;32m'
#YELLOW='\033[1;33m'
#BLUE='\033[0;34m'
#LBLUE='\033[1;34m'
#CYAN='\033[0;36m'
NC='\033[0m' # No Color => end color
#printf "I ${RED}love${NC} ${GREEN}Stack Overflow${NC}\n"

#---------------------------------------------------------------------------------#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION DEFINITIONS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#
#---------------------------------------------------------------------------------#

function print_start_time()
{
   echo -n "[$(date +%T)] "
}
#-----------------------------------------------------------------------------------------

function check_output()
{
    [ $DEBUG -eq 1 ] && msg " => working in $FUNCNAME ..." DEBUG NC
    outfile=$1
    pPID=$2

    if [ -s $outfile ]
    then
         msg " >>> wrote file $outfile ..." PROGR GREEN
	 return 0
    else
        echo
	msg " >>> ERROR! The expected output file $outfile was not produced, will exit now!" ERROR RED
        echo
	exit 9
	[ $DEBUG -eq 1 ] && msg "check_output running: kill -9 $pPID" DEBUG NC
	kill -9 $pPID
    fi
   [ $DEBUG -eq 1 ] && msg " <= exiting $FUNCNAME ..." DEBUG NC
}

#----------------------------------------------------------------------------------------- 


function check_dependencies()
{
    for programname in clustalo muscle
    do
       #if which $programname >/dev/null; then <== avoid which
       # see: http://stackoverflow.com/questions/592620/check-if-a-program-exists-from-a-bash-script

       bin=$(type -P $programname)
       if [ -z $bin ]; then
          echo
          printf "${RED}# ERROR: $programname not in place!${NC}\n"
          echo "# ... you will need to install \"$programname\" first or include it in \$PATH"
          echo "# ... exiting"
          exit 1
       fi
    done

    echo
    echo '# Run check_dependencies() ... looks good: all required binaries and perl scripts are in place.'
    echo
}
#----------------------------------------------------------------------------------------- 

function print_help()
{
   cat <<EOF
   $progname v.$VERSION usage:
   
   REQUIRED:
    -a <string> alignment algorithm [clustalo|muscle; default:$alignment_algorithm]
    -i <string> input fasta file name
    -h <FLAG> print this help
    -v <FLAG> print version
    -R <integer> RUNMODE
          1     standard msa
	  2     profile-profile alingnment
	  3     sequence to profile alignment
    
   OPTIONAL:
    
  
   NOTE1: XXX
   
   TODO: 

EOF

   check_dependencies
   
   exit 2  
}
#----------------------------------------------------------------------------------------- 


#------------------------------------#
#----------- GET OPTIONS ------------#
#------------------------------------#

input_fasta=
runmode=

alignment_algorithm=clustalo
DEBUG=0

# See bash cookbook 13.1 and 13.2
while getopts ':i:d:R:hD?:' OPTIONS
do
   case $OPTIONS in

   a)   alignment_algorithm=$OPTARG
        ;;
   i)   input_fasta=$OPTARG
        ;;
   h)   print_help
        ;;
   v)   echo "$progname v.$VERSION"
        ;;
   R)   runmode=$OPTARG
        ;;
   D)   DEBUG=1
        ;;
   \:)   printf "argument missing from -%s option\n" $OPTARG
   	 print_help
     	 exit 2 
     	 ;;
   \?)   echo "need the following args: "
   	 print_help
         exit 3
	 ;;
    *)   echo "An  unexpected parsing error occurred"
         echo
         print_help
	 exit 4
	 ;;	 
   esac >&2   # print the ERROR MESSAGES to STDERR
done

shift $((OPTIND - 1))

if [ -z "$input_fasta_extension" ]
then
       echo "# ERROR: no input fasta file extension defined!"
       print_help
       exit 1    
fi

if [ -z "$runmode" ]
then
       echo "# ERROR: no runmode defined!"
       print_help
       exit 1    
fi

if [ -z "$DEBUG" ]
then
     DEBUG=0 
fi


###>>> Exported variables !!!
#declare -x skip_seqs_gt=$skip_seqs_gt perl # export only2perl!!!  $ENV{skip_seqs_gt}

cat <<MSG 

 $progname vers. $VERSION run on $start_time with the following parameters: 
 wkdir=$wkdir | runmode=$runmode | DEBUG=$DEBUG 
 input_assembly=$input_fasta 
 alignment_algorithm=$aalignment_algorithm
 
 invocation: $progname ${args[*]}" PROGR YELLOW

MSG


#
# >>>> MAIN CODE <<<<
#

if [ "$runmode" -eq 1 ]
then
    echo "# XXX"

fi


