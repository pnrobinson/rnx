



#include "Nussinov.h"

#include <cstring>
#include <iostream>
#include <stack>
#include <algorithm>    // std::max




/*
 * The constructor allocates memory for the class variables 
 * rna_ (the original sequence)
 * mat_ (a two-dimensionary matrix for the Nussinov algorithm)
 * structure_ (a string of parentheses and dots representing the structure).
 */
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
  delete [] rna_;
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
      return 1;
  } else if (a=='U') {
    if (b=='G')
      return 1;
    if (b=='A')
      return 1;
  } else if (a=='A' && b=='U')
    return 1;
  else if (a=='C' && b=='G')
    return 1;
  
   return 0;
}



void
Nussinov::traceback2(int i, int j) {
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
    return traceback2(i+1,j);
  }
  if (cij== bonds(rna_[i],rna_[j]) + mat_[i+1][j-1]) { /* i is bonded with j */
    structure_[i]='(';
    structure_[j]=')';
    return traceback2(i+1,j-1);
  }
  for (int k=i+1;k<j;++k) {
    if (cij==mat_[i][k] + mat_[k+1][j]) {
      traceback2(i,k);
      traceback2(k+1,j);
    }
  }
}

void
Nussinov::traceback() {
  std::stack<std::pair<int,int> > stck;
  int len = get_len();
  memset(structure_, '.', len);
  structure_[len] = '\0';
  stck.push(std::make_pair(0, len-1));
  while (!stck.empty()) {
    int i = stck.top().first;
    int j = stck.top().second;
    stck.pop();
    std::cout << "i=" << i << " j=" << j << " bonds = " << bonds(rna_[i],rna_[j]) << std::endl;
    if (i>=j) {
      continue;
    } else if (mat_[i+1][j]==mat_[i][j]) {
      /* no match for base i */
      stck.push(std::make_pair(i+1,j));
    } else if (mat_[i][j-1] == mat_[i][j]) {
      /* no match for base j */
      stck.push(std::make_pair(i,j-1));
    } else if (mat_[i+1][j-1] + bonds(rna_[i],rna_[j]) == mat_[i][j]) {
      /* match base i<->j */
      std::cout << "match at i=" << i << " and j=" << j << std::endl;
      structure_[i]='(';
      structure_[j]=')';
      stck.push(std::make_pair(i+1,j-1));
    } else {
      for (int k=i+1;k<j;++k) {
	if (mat_[i][k] + mat_[k+1][j] == mat_[i][j]) {
	  /* match with two subsequences from (i..k)&(k+1..j) */
	  stck.push(std::make_pair(i,k));
	  stck.push(std::make_pair(k+1,j));
	  break;
	}
      }
    }
  }
}




void Nussinov::debugPrintMatrix() {
  std::cout << std::endl;
  //std::cout << "\t" << rna_[0];
  int len = get_len();
  for (unsigned int i=0;i<len;++i)
    std::cout << "\t" << rna_[i];
  std::cout << std::endl;

  std::cout << rna_[0];
  for (unsigned int j=0;j<len;++j) {
    std::cout << "\t" << mat_[0][j];
  }
  std::cout << std::endl;

  for (unsigned int i=1;i<len;++i) {
    std::cout << rna_[i];
    for (unsigned int j=0;j<i;++j)
      std::cout << "\t";
    for (unsigned int j=i;j<len;++j)
      std::cout << "\t" << mat_[i][j];
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

/*
 * The implementation of the Nussinov algorithm has been
 * adapted from Chapter 10.2 of Durbin et al., Biological
 * Sequence Analysis (1998) Cambridge University Press.
 * PSEUDOCODE   
 * for i = 2 to L do s(i,i-1) = 0; 
 * for i = 1 to L do s(i,i) = 0; 
 * for j = 2 to L do 
 *   for i = 1 to j-1 do 
 *     s(i,j) = max { 
 *           s(i+1,j), 
 *           s(i,j-1), 
 *           s(i+1,j-1) + d(i,j), 
 *           max i<k<j [s(i,k) + s(k+1,j)] 
 *          } 
 *
 */
const char * Nussinov::fold_rna() {
  unsigned int len = get_len();
  /* 1. Initialisation. (p. 270). Not needed because we initialise the
     entire matrix mat_ to 0. */
  /* 2. Recursion. Starting with all sequences of length 2 up to L */
  std::cout << rna_ << std::endl;
  for (unsigned int j=1;j<len;++j) {
    for (unsigned int i=0;i<j-1;++i) {
      int mx = std::max(mat_[i+1][j],mat_[i][j-1]); /* no bond for i or j*/
      std::cout << "i=" << i << " (" << rna_[i] << ") j="<<j<<" ("<< rna_[j] 
		<< "): bonds = " << bonds(rna_[i],rna_[j]) << std::endl;
      std::cout << "mx (no bond)=" << mx   << std::endl;
      std::cout << "bonds(rna_[i],rna_[j]) + mat_[i+1][j-1] =" << bonds(rna_[i],rna_[j]) + mat_[i+1][j-1]  << std::endl;
      mat_[i][j] = std::max(mx, bonds(rna_[i],rna_[j]) + mat_[i+1][j-1] ); /* bond i<->j */
      std::cout << "mat_[" << i << "][" << j<< "]=" << mat_[i][j] << std::endl;
      /* now check for pairs of subalignments */
      for (int k=i+1;k<j;++k) {
	/* now check for pairs of subalignments */
	mat_[i][j] = std::max(mat_[i][j],
			      mat_[i][k] + mat_[k+1][j]);
      }
    }
  }
  // When we get here the alignment is finished and we need to do the traceback
  traceback();
  debugPrintMatrix();
  return structure_;
}

