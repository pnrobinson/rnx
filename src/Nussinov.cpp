



#include "Nussinov.h"

#include <cstring>
#include <algorithm>    // std::max





Nussinov::Nussinov(const char * rna) {
  rna_ =  new char[strlen(rna) + 1]; 
  strcpy(rna_, rna);
}

Nussinov::~Nussinov() {
  delete rna_;
  rna_=NULL;
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
  else if (a=='G' && b=='C')
    return 3;
  else
    return 0;
}



void traceback(int i, int j, char *p) {


}


char * Nussinov::fold_rna() const {
 
  unsigned int len = get_len();
  int** ary = new int*[len];
  for(int i = 0; i < len; ++i) {
    ary[i] = new int[len];
    memset(&(ary[i]),0,len*sizeof(int));
  }
  // p will store the paren-dot string representing the structure.
  char *p = new char[len+1];
  memset(p, '\0', len+1 );
  // No hydrogen bonds are possible in less than 5 base pairs
  // start loop at fifth position (i=4)
  for (int z=4;z<len;++z) {
    for (int i=0;i<(len-z+1);++i) {
      int j=i+z-1;
      ary[i][j]= std::max(ary[i+1][j] /* no bond */,
			  bonds(rna_[i],rna_[j]) +
			  ary[i+1][j-1] /* bond i<->j */);
      for (int k=i+1;k<j;++k) {
	/* now check for pairs of subalignments */
	ary[i][j] = std::max(ary[i][j],
			     ary[i][k] + ary[k+1][j]);
      }
    }
  }
  // When we get here the alignment is finished and we need to do the traceback

  // clean up
  for(int i = 0; i < len; ++i) {
    delete [] ary[i];
  }
  delete [] ary;
  return p;
}

