/*
  NAME SVDMAN.CC

  AUTHORS

     Michael Wall, Los Alamos National Laboratory/mewall@lanl.gov
     Patricia Dyck, Los Alamos National Laboratory/tdyck@lanl.gov

  SVDMAN (svdman.cc, svdtool.cc, svd.h) was
  developed at Los Alamos National Laboratory under the auspices of
  the Department of Energy under contract to the University of
  California.

  Please see html documentation in ../doc/svdman.html

  Copyright (c) 2000,2001 The Regents of the University of California.
  All rights reserved.

  This software and ancillary information (herein called "SOFTWARE") called
  SVDMAN 0.9b is made available under the terms described here. The SOFTWARE 
  has been approved for release with associated LA-CC number 01-12.

  Unless otherwise indicated, this SOFTWARE has been authored by an employee 
  or employees of the University of California, operator of the Los Alamos 
  National Laboratory under Contract No. W-7405-ENG-36 with the US Dept 
  of Energy. The US Government has rights to use, reproduce, and distribute 
  this SOFTWARE. The public may copy, distribute, prepare derivative works 
  and publicly display this SOFTWARE without charge, provided that this 
  Notice and any statement of authorship are reproduced on all copies. 
  Neither the Government nor the University makes any warranty, express or 
  implied, or assumes any liability or responsibility for the use of this 
  SOFTWARE.

  If SOFTWARE is modified to produce derivative works, such modified 
  SOFTWARE should be clearly marked, so as not to confuse it with the 
  version available from LANL.

*/

#include "svd.h"

int main(int argc, char *argv[]){

  ios::sync_with_stdio();

  std::string ifname, out_base;

  double thresh=3.0;

  // Flags to direct actions (0 means "no"):

  char 
    have_ifname = 0,
    have_out_base = 0, 
    have_minus = 0, 
    have_submean = 0,
    have_writesvd = 0,
    have_evaluate = 0;

  /*
   * Parse command line arguments: */

  int inpct=0;
  char *cmd;

  --argc; 
  while (argc>0) {
    if (argv[argc][0] != '-') {
      ++inpct;
    }
    else {
      cmd = argv[argc]+1;
      /*
       * Input file name:
       */
      if (!strcmp(cmd,"i")) {
	if (inpct == 1) {
	  ifname = argv[argc+1];
	  have_ifname = 1;
	}
	else {
	  cerr << "SVDMAN_ERR: -i takes 1 argument only.\n";
	  exit(1);
	}
	inpct=0;
      }
      /*
       * Output base name:
       */
      if (!strcmp(cmd,"o")) {
	if (inpct == 1) {
	  out_base = argv[argc+1];
	  have_out_base = 1;
	}
	else {
	  cerr << "SVDMAN_ERR: -o takes 1 argument only.\n";
	  exit(1);
	}
	inpct=0;
      }
      /*
       * Threshold for inclusion in groups:
       */
      if (!strcmp(cmd,"thresh")) {
	if (inpct == 1) {
	  thresh = atof(argv[argc+1]);
	}
	else {
	  cerr << "SVDMAN_ERR: -t takes 1 argument only.\n" ;
	  exit(1);
	}
	inpct=0;
      }
      /*
       * Flag to reverse signs of VT and U:
       */
      if (!strcmp(cmd,"minus")) {
	if (inpct == 0) {
	  have_minus = 1;
	}
	else {
	  cerr << "SVDMAN_ERR: -minus takes 0 arguments only.\n";
	  exit(1);
	}
	inpct=0;
      }
      /*
       * Flag to subtract the mean value from each row of A:
       */
      if (!strcmp(cmd,"submean")) {
	if (inpct == 0) {
	  have_submean = 1;
	}
	else {
	  cerr << "SVDMAN_ERR: -submean takes 0 arguments only.\n";
	  exit(1);
	}
	inpct=0;
      }
      /*
       * Flag to output the SVD in four files: U,S,VT,B
       */
      if (!strcmp(cmd,"writesvd")) {
	if (inpct == 0) {
	  have_writesvd = 1;
	}
	else {
	  cerr << "SVDMAN_ERR: -writesvd takes 0 arguments only.\n";
	  exit(1);
	}
	inpct=0;
      }
      /*
       * Flag to calculate self-consistency
       */
      if (!strcmp(cmd,"evaluate")) {
	if (inpct == 0) {
	  have_evaluate = 1;
	}
	else {
	  cerr << "SVDMAN_ERR: -evaluate takes 0 arguments only.\n";
	  exit(1);
	}
	inpct=0;
      }
    }
    --argc;
  }

  if (!have_ifname) {
    cerr <<
	    "Usage:" << endl <<

        "   Execution:" << endl <<

        "       svdman -i [ifname] { -options {[args]} }" << endl <<

        "   Commands:" << endl <<

        "       Required:" << endl << 

        "           -i [std::string ifname]" << endl << 

        "              ifname: input file name" << endl << 

        "       Options:" << endl <<  

        "           -o [std::string ofbase]" << endl << 

        "              ofbase: output filename base" << endl << 

        "           -thresh [double threshold]" << endl << 

        "              threshold: A scale for the threshold value beneath which an ORF does not belong to an eigenvector" << endl << 

        "           -writesvd: flag to write SVD components U, S, VT and B to files." << endl <<

        "           -minus: flag to invert the sign of the matrix elements for U, VT and B" << endl <<

        "           -submean: flag to subtract the row mean value from each row" << endl <<

        "           -evaluate: flag to check for self-consistency using an interpolation procedure" << endl << endl;
    exit(1);
  }

  // Read SVD data from file into d:

  Data d(ifname);

#ifdef DEBUG
  cerr << "Data read from file " << ifname << endl;
#endif

  // Subtract the mean of the row from each data value:

  if (have_submean) {
    d.SubMeanRow();
  }

  // Calculate the SVD:

  d.SVD();

#ifdef DEBUG
  cerr << "SVD calculated.\n";
#endif

  // Invert the sign of the matrix elements for U, VT if desired:

  if (have_minus) {
    d.MinusSVD();
  }

  // Output groupings if an output base is specified:

  if (have_out_base) {
    d.WriteGroups(out_base, thresh);
  }

  //  Output the SVD results if desired:
  
  if (have_writesvd) {
    d.WriteSVD(out_base);
  }

  // Perform self-consistency validation if desired:

  if (have_evaluate) {
  
    std::vector<double> correlations;

#ifdef DEBUG
    cerr << "Validating...\n";
#endif

    d.EvaluateSVD(correlations);

#ifdef DEBUG
    cerr << "   evaluated.\n";
#endif

    if (have_out_base) {
      d.WriteCorr(out_base, correlations);
    }

  } // if have_evaluate 

} //main

