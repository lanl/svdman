/*
  NAME SVDTOOL.CC

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

int WriteRawData(std::string ofname, std::vector<double> dat, int numcols, std::vector<std::string> row_labels, std::vector<std::string> col_labels) {
  ofstream data_file;

  data_file.open(ofname.c_str());
  if (data_file.bad()) {
    cerr << "WriteRawData_ERR: Can't open " << ofname << ".\n";
    exit(1);
  }
  
  for (int i=0;i<col_labels.size();i++) {
    data_file << '\t' << col_labels[i];
  }
  data_file << "\n";

  for (int i=0;i<row_labels.size();i++) {
    data_file << row_labels[i];
    for (int j=0;j<col_labels.size(); j++) {
      data_file << "\t" << dat[i*numcols + j];
    }
    data_file << "\n";
  }
  data_file.close();
}

void splint(double xa[],double ya[],double y2a[],int n,double x, double& y)
{
	int klo,khi,k;
	double h,b,a;
	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k-1] > x) khi=k;
		else klo=k;
	}
	khi--;
	klo--;
	h=xa[khi]-xa[klo];
	if (h == 0)
	cout << "h = 0" << endl;
	/*if (h == 0.0) nrerror("Bad xa input to routine splint");*/
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;


	y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

void spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
{
	int i, k ;
	double p,qn,sig,un;

	double *u=(double *)calloc(n,sizeof(double));

	if (yp1 > 0.99e30)
		y2[0]=u[0]=0.0;
	else {
		y2[0] = -0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
	
	for (i=1;i<=n-2;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}

	
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	
	
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k=n-2;k>=0;k--){
		y2[k]=y2[k]*y2[k+1]+u[k];
	}
	free((double *)u);
}

int Max(int m, int n)
{
  int rval;

  if (m>n) {rval=m;} else {rval=n;}

  return(rval);
}

int Min(int m, int n)
{
  int rval;

  if (m<n) {rval=m;} else {rval=n;}

  return(rval);
}

int SplitLine(std::string strline, int start, int nwords, char ch, std::vector<std::string> &words) {

  unsigned int pos=0;
  int oldpos=0, i=0, nfound = 0;

  int stop=start+nwords;

  while (((pos=strline.find(ch,oldpos)) < strline.length()) && 
	 ((i < stop ) || (nwords == 0))) {
    if (i>=start) {
      words.push_back(strline.substr(oldpos,pos-oldpos));
      nfound++;
    }
    oldpos = pos+1;
    i++;
  }
  if ((i < stop) || (nwords == 0)) {
    words.push_back(strline.substr(oldpos,strline.length()-oldpos));
    nfound++;
  }

  return(nfound);
}

Data::Data(std::string ifname)
{
  ifstream data_file;
  std::string input_line;

  int numsplit;

  std::string valstr;

  std::vector<std::string> words;

  data_file.open(ifname.c_str());
  if (data_file.bad()) {
    cerr << "ReadDataERR: Can't open " << ifname << ".\n";
    exit(1);
  }

  // Get column labels:

  getline(data_file,input_line);
  SplitLine(input_line,1,0,'\t',col_labels);
  numcols = col_labels.size();
  int rownum=0;
  while (getline(data_file,input_line)) {
    SplitLine(input_line,0,1,'\t',row_labels);
    words.clear();
    numsplit = SplitLine(input_line,1,numcols,'\t',words);
    if (numsplit<numcols) {
      cerr << "ReadDataERR: Number of fields is incorrect for row " << rownum 
	   << ".\n";
      exit(1);
    }
    for (int i=0;i<words.size();i++) {
      data.push_back(atof(words[i].c_str()));
    }
    rownum++;
  }
  numcols = col_labels.size();
  numrows = row_labels.size();
  data_file.close();
}

Data::Data(int index, Data& d)
{
  numrows = d.numrows;
  numcols = d.numcols -1;
  for (int i=0;i < d.numrows; i++) {
    row_labels.push_back(d.row_labels[i]);
  }
  //copy col_labels
  for (int i = 0; i < d.numcols; i++) {
    if (i != index) {
      col_labels.push_back(d.col_labels[i]);
    }//for
  }
  for (int j = 0; j < d.numrows; j++) {
    for (int i = 0; i < d.numcols; i++) {
      if (i != index) {
	data.push_back(d.data[j*d.numcols+i]);
      } //if
    } //for		 
  } //for
#ifdef DEBUG1
  cerr << " Size of reduced data array = " << data.size() << "\n";
#endif
}

Data::Data(Data& d)
{
  numrows = d.numrows;
  numcols = d.numcols;
  row_labels.clear();
  for (int i=0;i < d.numrows; i++) {
    row_labels.push_back(d.row_labels[i]);
  }
  //copy col_labels
  col_labels.clear();
  for (int i = 0; i < d.numcols; i++) {
    col_labels.push_back(d.col_labels[i]);
  }
  data.clear();
  for (int j = 0; j < d.numrows; j++) {
    for (int i = 0; i < d.numcols; i++) {
      data.push_back(d.data[j*d.numcols+i]);
    } //for		 
  } //for
}

int Data::NumCols() {return numcols;};

int Data::NumRows() {return numrows;};

int Data::WriteCorr(std::string out_base, std::vector<double> correlations) {

  std::string ofname = out_base + "_evaluate.dat";
  ofstream data_file;
  
  data_file.open(ofname.c_str());
  if (data_file.bad()) {
    cerr << "SVDMAN_ERR: Can't open " << ofname << ".\n";
    exit(1);
  }
  
  // Output correlations
  data_file << "\tCorrelation\n";
  for (int i = 0 ; i < correlations.size(); i++){
    data_file << col_labels[i+1] << "\t" << correlations[i] << endl;
  }
  
  data_file.close();

}  
  
int Data::WriteGroups(std::string out_base, double thresh) {

  ofstream data_file;
  double t = thresh/sqrt((double)row_labels.size());

  std::string ofname = out_base + "_groups.dat";

  data_file.open(ofname.c_str());
  if (data_file.bad()) {
    cerr << "WriteRawData_ERR: Can't open " << ofname << ".\n";
    exit(1);
  }
  
  for (int i=0;i<mode_labels.size();i++) {
    for (int j=0;j<row_labels.size(); j++) {
      if (fabs(VT[j*numcols + i]) > t) {
	data_file << mode_labels[i] << "\t" << row_labels[j] << "\t";
	if (VT[j*numcols + i]<0) {
	  data_file << "-" << "\n";
	} else {
	  data_file << "+" << "\n";
	}
      }
    }
  }
  data_file.close();
}

int Data::SubMeanRow()
{
  for (int i=0;i<numrows;i++) {
    double mean=0;
    for (int j=0; j<numcols;j++) {
      mean += data[i*numcols + j];
    }
    mean /= (double)numcols;
    for (int j=0; j<numcols;j++) {
      data[i*numcols + j] -= mean;
    }
  }
  return 0;
}

void Data::SVD(){

  //set stuff needed for SVD
  
  int Alen = numrows*numcols;

  double *A = (double *)calloc(Alen,sizeof(double));

  if (!A) {
    cerr << "DATA::SVD_ERR: Couldn't allocate array A.\n";
    exit(1);
  }

  for (int i=0;i< data.size();i++) {
    A[i] = data[i];
  }

  // Allocate arrays Utmp, VTtmp, Stmp and W:

  int Ulen = numcols*numcols;

  double *Utmp = (double *)calloc(Ulen,sizeof(double));

  if (!Utmp) {
    cerr << "DATA::SVD_ERR: Couldn't allocate array Utmp.\n";
    exit(1);
  }

  // Don't actually use the array VT.  
  // Instead return VT in A.  
  // Create dummy pointer here.

  int VTlen = 1;

  double *VTtmp = (double *)calloc(VTlen,sizeof(double));

  if (!VTtmp) {
    cerr << "DATA::SVD_ERR: Couldn't allocate array VTtmp.\n";
    exit(1);
  }

  int Slen = numcols;

  double *Stmp = (double *)calloc(Slen,sizeof(double));

  if (!Stmp) {
    cerr << "DATA::SVD_ERR: Couldn't allocate array Stmp.\n";
    exit(1);
  }

  int Wlen = 10*Max(3*Min(numrows,numcols)+Max(numrows,numcols),5*Min(numrows,numcols)-4);

  double *W = (double *)calloc(Wlen,sizeof(double));

  if (!W) {
    cerr << "DATA::SVD_ERR: Couldn't allocate array W.\n";
    exit(1);
  }

  // SVD subroutine arguments (both ask for concise calculation):

  char jobu = 'S';  // return U in Utmp
  char jobvt = 'O'; // return VT in overwritten A

  // SVD return value:

  int info;

  // Calculate the SVD:

  dgesvd_(&jobu,&jobvt,&numcols,&numrows,A,&Slen,Stmp,Utmp,&numcols,VTtmp,&numrows,W,&Wlen,&info);

#ifdef DEBUG1
  cerr << "Made it past SVD.\n";
#endif

  // Copy results into working arrays:

  U.clear();

  for (int i=0;i < Ulen;i++) {
    U.push_back(Utmp[i]);
  }
  free((double *)Utmp);

  VTlen = numcols*numrows;

  VT.clear();

  for (int i=0;i < numrows;i++) {
    for (int j=0;j < numcols;j++) {
      VT.push_back(A[i*numcols+j]);
    }
  }
  free((double *)VTtmp);
  free((double *)A);

#ifdef DEBUG1
  cerr << "Made it past transfer to working arrays.\n";
#endif 

  S.clear();

  for (int i=0;i < Slen;i++) {
    S.push_back(Stmp[i]);
  }
  free((double *)Stmp);

  // Create array B:

  B.clear();

  for (int i  = 0 ; i < numcols; i ++){
    for (int j = i*numcols ; j < (i*numcols)+numcols ; j++){
      B.push_back(U[j] * S[i]);
    }
  }	

  mode_labels.clear();
  for (int i=0;i<numcols;i++) {
    char s[20];
    sprintf(s,"%03d",i);
    mode_labels.push_back(s);
  }

}

void Data::MinusSVD(){

  for (int i=0;i<VT.size();i++) {
    VT[i] = -VT[i];
  } 

  for (int i=0;i<U.size();i++) {
    U[i] = -U[i];
    B[i] = -B[i];
  }

}

void Data::WriteSVD(std::string out_base){

  // ***NOTE: VT and U are swapped here in order to conform to a standard where
  //          the genes are indexed by rows, assays by columns.

  std::string fname = out_base + "_VT.dat";
    
  WriteRawData(fname,U,col_labels.size(),mode_labels,col_labels);
  
  fname = out_base + "_U.dat";
  
  WriteRawData(fname,VT,mode_labels.size(),row_labels,mode_labels);
  
  std::vector<std::string> s_label;
  s_label.push_back("Value");
  
  fname = out_base + "_S.dat";
  
#ifdef DEBUG 
  cerr << "Writing singular values S.\n";
#endif
  WriteRawData(fname,S,mode_labels.size(),s_label,mode_labels);
#ifdef DEBUG 
  cerr << "  Wrote singular values S.\n";
#endif
  
  fname = out_base + "_B.dat";
  
  WriteRawData(fname,B,col_labels.size(),mode_labels,col_labels);
  
//  cout << "Number of targets (length of singular std::vectors): " << numrows << "\n"
 //      << "Number of assays (number of singular std::vectors): " << numcols << "\n";
  
}	//writesvd


int Data::JacknifeData(double xval, std::vector<double> &data_col)
{
    std::vector <double> Bcol;
    InterpolateSVD(xval,Bcol);
    CalcColData(Bcol, data_col);
    return 0;
}

int Data::EvaluateSVD(std::vector<double> &correlations)
{

  for (int i = 1 ; i < numcols-1 ; i++){
    
#ifdef DEBUG1
    cerr << i << "...";
#endif
    
    //create new data struct and erase column i
    
    Data new_d(i, *this);
    
#ifdef DEBUG1
    cerr << "Creating subarray.\n";
#endif
    
#ifdef DEBUG1
    cerr << "  Created subarray.\n";
#endif
    
#ifdef DEBUG1
    cerr << "  Calculating SVD of subarray.\n";
#endif
    
    new_d.SVD();
    
#ifdef DEBUG1
    cerr << "  Calculated SVD of subarray.\n";
    new_d.writeSVD("test",thresh);
#endif
    
    std::vector<double> data_col;

    new_d.JacknifeData((double)atof(col_labels[i].c_str()),data_col);

    
#ifdef DEBUG1
    cerr << "  Interpolated values from subarray.\n";
#endif
    
    //make new std::vector to store calculated A values
    
#ifdef DEBUG1
    cerr << "data_col.size() = " << data_col.size() << "\n";
#endif
    
    CorrColData(i, data_col, correlations);
    
  }	// end of for loop
  
  return 0;

}
  
int Data::CalcColData(std::vector<double>& calcValues, std::vector <double> &data_col)
{
  double sum;
// in each row
  for (int i = 0; i < numrows;i++) {
    //initialize sum to zero
    sum = 0;
    //go through the row 
    for (int j = 0 ; j < numcols; j++){
      sum = sum + (VT[i*numcols + j]*calcValues[j]);
    } 
    data_col.push_back(sum);
  }
  return 0;
}                                   

int Data::InterpolateSVD(double xval, std::vector<double> &newValues)
{
  double *column_values = (double *)malloc((numcols+1)*sizeof(double));
  double *row_values = (double *)malloc((numcols+1)*sizeof(double));
  double *y2 = (double *)malloc((numcols+1)*sizeof(double));
 
  
  for (int i=0;i<numcols;i++) {
    column_values[i] = atof(col_labels[i].c_str());
  }

  // go through each line of the data, writing data to an array 
  // interpolating, and storing the new value in std::vector
  
  for (int i = 0; i < numcols ; i++){
    for (int j = 0 ; j < numcols ; j++){
      row_values[j] = B[i*numcols + j];	
    }	//for

    double value;

    spline (column_values, row_values, numcols, 1.0e30, 1.0e30, y2);
    splint (column_values, row_values, y2, numcols, xval, value);


#ifdef DEBUG1
    cerr << "std::vector " << i << " value = " << value << "\n";
#endif

    newValues.push_back(value);
  }	//for

  free((double *)column_values);
  free((double *)row_values);
  free((double *)y2);

  return 0;
}

void Data::CorrColData(int index, std::vector<double>& data_col, std::vector<double>& correlations)
{
  double Xtot, Ytot, Xavg, Yavg, numerator, Xdenominator, Ydenominator;
  Xtot = 0;
  Ytot = 0;
  numerator = 0;
  Xdenominator = 0;
  Ydenominator = 0;
  
  // get X avg
  for (int j = 0 ; j < numrows ; j++){
    Xtot = Xtot + data[j*numcols + index];
  }
  Xavg = Xtot / (double)numrows;

  
  // get Y avg
  for (int i = 0 ; i < numrows ; i++){
    Ytot = Ytot + data_col[i];
  }
  Yavg = Ytot / (double)numrows;
  
#ifdef DEBUG1
  cerr << "Xavg = " << Xavg << "; Yavg = " << Yavg << "\n";
#endif

  for (int j = 0 ; j < numrows ; j++){
    numerator = numerator + ((data[j*numcols+index] - Xavg) * (data_col[j] - Yavg));
    Xdenominator = Xdenominator + ((data[j*numcols+index] - Xavg)*
				   (data[j*numcols+index] - Xavg));
    Ydenominator = Ydenominator + ((data_col[j] - Yavg)*
				   (data_col[j] - Yavg));
  }
  
  double denominator = sqrt( Xdenominator * Ydenominator) ;
  
  double number = numerator / denominator;
  
  correlations.push_back(number);
}
