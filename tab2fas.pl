#!/usr/bin/perl 

# Pablo Vinuesa
# Centro de Ciencias Genomicas, UNAM, Mexico
# Jul 26th, 2010

# tab2fas, restores a fasta file from a fastab-formated file
# translated from the onliner: 
# perl -pe 'if(/>/){s/\t/\n/}' infile.fastab > outfile.fas

use File::Basename;

$progname = basename($0); #tab2fas
$version = 1.0;

die "\n# $progname v.$version needs a fastab (tabular fasta) file name provided as single argument\n"
  unless @ARGV == 1;
while(<>)
{
   s/^\n$//; # get rid of the 1st empty line produced by fas2tab.pl
   if(/>/)
   {
       s/\t/\n/;
   }

}
continue
{
   print
}
