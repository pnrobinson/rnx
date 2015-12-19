#include "Nussinov.h"

#include <cstring>
#include <iostream>
#include <stack>
#include <algorithm>    // std::max




/**
 * The constructor allocates memory for the class variables 
 * rna_ (the original sequence)
 * mat_ (a two-dimensionary matrix for the Nussinov algorithm)
 * structure_ (a string of parentheses and dots representing the structure).
 * @param rna a string (upper case) of RNA nucleotides 
 */
Nussinov::Nussinov(const char * rna): h_(3) {
  init(rna);
}

/**
 * This constructor allows client code to adjust the minimum
 * distance between two baired based (\code{h_}).
 */
Nussinov::Nussinov(const char * rna, unsigned int h): h_(h) {
  init(rna);
}

/**
 * Allocate memory for the constructors
 * @param rna a string (upper case) of RNA nucleotides 
 */
void Nussinov::init(const char * rna) {
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

/** @return length of the RNA sequence being analysed. */
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



/**
 * Prints the DP matrix to the shell, and can
 * be used for debugging with relatively short
 * RNA sequences.
 */
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

/**
 * Return a parenthesis-dot representation of the RNA secondary structure. 
 * Note that clients should make a copy of this string if they need to alter it
 * The implementation of the Nussinov algorithm has been
 * adapted from Chapter 10.2 of Durbin et al., Biological
 * Sequence Analysis (1998) Cambridge University Press.
 * Note that the following pseudocode uses 1-based numbering. The for loops
 * iterate through each diagonal above the main diagonal one at a time, ensuring
 * that all information needed for the Dynamic Programming step is available when needed.
 * \code
 * for i = 2 to L
 *   s[i,i-1] = 0; // initialize subdiagonal 
 * for i = 1 to L
 * s[i,i] = 0;     // initialize main diagonal 
 * for s = 2 to L
 *   i = 0;  
 *   for j = s to L
 *     i=i+1
 *     s[i,j) = max { 
 *           s(i+1,j), 
 *           s(i,j-1), 
 *           s(i+1,j-1) + d(i,j), 
 *           max i<k<j [s(i,k) + s(k+1,j)] 
 *          }
 * \endcode
 *
 */
const char * Nussinov::fold_rna() {
  unsigned int len = get_len();
  /* 1. Initialisation. (p. 270). Not needed because we initialise the
     entire matrix mat_ to 0. */
  /* 2. Recursion. Starting with all sequences of length 2 up to L */
  /* Note we start at the subdiagonal after h_ bounds, this prevents
     bairs from being taken that are too close together. */
  for (unsigned int s=1+h_;s<len;++s) {
    unsigned int i=0;
    for (unsigned int j=s;j<len;++j,++i) {
   
      int mx = std::max(mat_[i+1][j],mat_[i][j-1]); /* no bond for i or j*/
      mat_[i][j] = std::max(mx, bonds(rna_[i],rna_[j]) + mat_[i+1][j-1] ); /* bond i<->j */
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
  //debugPrintMatrix();
  return structure_;
}


