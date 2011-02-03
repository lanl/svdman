/*
  NAME SVD.H

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

#include "stdlib.h"
#include "f2c.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include<string>
#include<vector>
#include "iostream.h"
#include "fstream.h"

class Data{

 private:

  // Variables:

  std::vector<double> data; // Original data array with gene expression levels
  int numrows; // Number of rows in data
  int numcols; // Number of columns in data
  std::vector<std::string> row_labels; // Row labels
  std::vector<std::string> col_labels; // Column labels
  std::vector<std::string> mode_labels; // Mode labels for SVD (indexes)
  std::vector<double> B; // Matrix U*S
  std::vector<double> VT; // Matrix of right singular vectors (indexed by column)
  std::vector<double> S; // List of singular values
  std::vector<double> U; // Matrix of left singular vectors (indexed by row)

 public:

  // Copy constructors:

  Data (std::string ifname); // Read from a file
  Data (int index, Data& original); // Derive from existing w/o column "index"
  Data (Data& original); // Copy from existing

  // Status reports:

  int NumRows(); // Return the number of rows
  int NumCols(); // Retun the number of columns

  // Input/Output:

  int WriteGroups(std::string out_base, double thresh); // Write SVD-based groupings
  void WriteSVD(std::string out_base); // Write the SVD
  int WriteCorr(std::string out_base, 
		std::vector<double> correlations); // Write the correlations

  // Transformations:

  int SubMeanRow(); // Subtract from each element the mean row value
  void MinusSVD(); // Invert the sign of the elements of U and VT

  // Calculations:

  void SVD(); // SVD calculateion
  int EvaluateSVD(std::vector<double> &correlation); // Correlations from valid.
  int JacknifeData(double xval, std::vector<double> &A_calc);
  void CorrColData(int index, std::vector<double>& data_col, std::vector<double>& correlations);
  int InterpolateSVD(double xval, std::vector<double> &Bcol);
  int CalcColData(std::vector<double>& calcValues, std::vector<double> &data_col);
};

extern "C" int dgesvd_(char *jobu, char *jobvt, int *m, int *n, 
 	    double *a, int *lda, double *s, double *u, int *
	    ldu, double *vt, int *ldvt, double *work, int *lwork, 
  			     int *info);
void splint(double xa[],double ya[],double y2a[],int n,double x, double& y);
void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
int Max(int m, int n);
int Min(int m, int n);
int SplitLine(std::string strline, int start, int nwords, char ch, std::vector<std::string> &words);
int WriteRawData(std::string ofname, std::vector<double> dat, int numcols, std::vector<std::string> row_labels, std::vector<std::string> col_labels);
