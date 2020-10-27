#!/usr/bin/perl 

# Pablo Vinuesa
# Centro de Ciencias Genomicas, UNAM, Mexico
# Jul 26th, 2010

# fas2tab, translated from the onliner: 
# perl -pe 'unless(/^>/){s/\n//g}; if(/>/){s/\n/\t/g}; s/>/\n>/' infile > outfile.fastab

use File::Basename;

$progname = basename($0); #fas2tab
$version = 1.0;

die "\n# $progname v.$version needs a fasta file name provided as single argument\n"
  unless @ARGV == 1;
while(<>)
{
   unless(/^>/)
   {
       s/\n//g
   }
   if(/>/)
   {
       s/\n/\t/g
   }
   s/>/\n>/;
}
continue
{
   print
}
