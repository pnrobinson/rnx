



#include "Nussinov.h"

#include <cstring>
#include <iostream>
#include <algorithm>    // std::max





Nussinov::Nussinov(const char * rna) {
  unsigned int len = strlen(rna);
  rna_ =  new char[len + 1];
  strcpy(rna_, rna);
  mat_ = new int*[len];
  for(int i = 0; i < len; ++i) {
    mat_[i] = new int[len];
    memset(mat_[i],0,len*sizeof(int));
  }
  structure_ = new char[len+1];
  memset(structure_, '\0', len+1 );
}

Nussinov::~Nussinov() {
  unsigned int len = get_len();
  delete rna_;
  rna_=NULL;
  for(int i = 0; i < len; ++i) {
    delete [] mat_[i];
  }
  delete [] mat_;
  delete [] structure_;
}

unsigned int Nussinov::get_len() const {
  return strlen(rna_);
}


/**
 * @return number of hydrogen bonds between RNA nucleotides
 * a and b if  a and b are paired in the RNA secondary structure.
 */
int bonds(char a, char b) {
  if (a=='G'){
    if (b=='U')
      return 1;
    if (b=='C')
      return 3;
  } else if (a=='U') {
    if (b=='G')
      return 1;
    if (b=='A')
      return 2;
  } else if (a=='A' && b=='U')
    return 2;
  else if (a=='C' && b=='G')
    return 3;
  
   return 0;
}



void
Nussinov::traceback(int i, int j) {
  unsigned int len = get_len();
  if (i>=j || i>=len)
    return;
  int cij=mat_[i][j];
  if (cij==0) {
    for (int k=i;k<=j;++k)
      structure_[k]='.';
  }
  if (cij==mat_[i+1][j]) { /* position i is unbonded */
    structure_[i]='.';
    return traceback(i+1,j);
  }
  if (cij== bonds(rna_[i],rna_[j]) + mat_[i+1][j-1]) { /* i is bonded with j */
    structure_[i]='(';
    structure_[j]=')';
    return traceback(i+1,j-1);
  }
  for (int k=i+1;k<j;++k) {
    if (cij==mat_[i][k] + mat_[k+1][j]) {
      traceback(i,k);
      traceback(k+1,j);
    }
  }
}


char * Nussinov::fold_rna() {
 
  unsigned int len = get_len();
  // No hydrogen bonds are possible in less than 5 base pairs
  // start loop at fifth position (i=4)
  for (int z=4;z<len;++z) {
    for (int i=0;i<(len-z);++i) {
      int j=i+z;
      mat_[i][j]= std::max(mat_[i+1][j] /* no bond */,
			  bonds(rna_[i],rna_[j]) +
			  mat_[i+1][j-1] /* bond i<->j */);
      for (int k=i+1;k<j;++k) {
	/* now check for pairs of subalignments */
	mat_[i][j] = std::max(mat_[i][j],
			     mat_[i][k] + mat_[k+1][j]);
      }
    }
  }
  // When we get here the alignment is finished and we need to do the traceback
  traceback(0,len-1);
 
  return structure_;
}

