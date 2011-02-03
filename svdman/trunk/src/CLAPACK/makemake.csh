#!/bin/csh -f
#
# MAKEMAKE - Make the master Makefile for a c-code package
# Author: Mike Wall
# Date: 11/21/95
# Description:  
#               When in the ./c/src directory, this creates a makefile which
#            builds the executables.
# Usage:
#               makemake <makefile>
#
# Save the old one:
#
cp $1 $1.sav
#
# Get the header:
#
cat $1.head > $1
#
# Primary list:
#
echo "all: \" >> $1
ls -1 *.c | gawk '{print "\t"substr($0,0,length()-2)".o \\"}' >> $1
echo "\tlibclasvd.a\n" >> $1
#
# Build library functions:
#
echo "#\n# Library functions:\n#" >> $1
ls -1 *.c | gawk '{print substr($0,0,length()-2)".o:"$1"\n\t$(CC) -c -o "substr($0,0,length()-2)".o $(CFLAGS) "$0}' >> $1
#
# Build libraries:
#
echo "#\n# Libraries:\n#\nlibclasvd.a: \" >> $1
ls -1 *.c | awk '{print "\t"substr($0,0,length()-2)".o \\"}' >> $1
echo -n '\t;$(AR) $(ARFLAGS) ../libclasvd.a ' >> $1
ls -1 *.c | awk '{ORS=" ";print substr($0,0,length()-2)".o"}' >> $1
echo '\n' >> $1

