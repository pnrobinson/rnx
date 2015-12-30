#ifndef RNA_STRUCTURE_H
#define RNA_STRUCTURE_H



#include <string>


#if !defined(DEFINES_H)
#define DEFINES_H
#define maxfil 100    //maximum length of file names
#define infinity 9999999  //an arbitrary value given to infinity
#define maxtloop 100 //maximum tetraloops allowed (info read from tloop)
#define maxstructures 1010 //maximum number of structures in ct file

#define ctheaderlength 125 //max length of string containing info on sequence
#define ga_bonus -10 //bonus value for "almost coaxial stacking" in efn2
#define amax 400 //this is a maximum line length for void linout (below)
#define col 80  //this is the number of columns in an output file
#define numlen 8  //maximum digits in a number
#define maxforce 600 //maximum number of bases that can be forced single
#define maxgu 5 //maximum number of u's in gu pair
#endif

/**
 * \class RNAStructure
 *
 * \ingroup Folding
 * 
 * This class represents a single RNA Structure
 * that was inferred by an algorithm or input from
 * a file such as a CT file.
 * This class et up to hold many possible structures of the same sequence
 * @author Peter Robinson
 * @version 0.0.3 Dec 30, 2015
 */
class RNAStructure {
  /** (Maximum) length of the header line of a CT file. */
  static const int s_ctheaderlength = 125;
  /** maximum number of structures in ct file */
  static const int s_maxstructures = 1010;
  /** maximum number of bases in a structure */
  static const int s_maxbases = 10000;
  /** number of bases in sequence */
  int numofbases_;
  /** number of alternative structures of the sequence
      that is held by structure */
  int numofstructures_;
  int pair[maxforce][2];
  int npair_;
  /**
   *  numseq[i] = a numeric that stands for the base in the ith position
    of the sequence,
    A = 1
    C = 2
    G = 3
    U = 4
  */
  int *numseq_;
  /** array stores the historical numbering of a sequence */
  int *hnumber_;
  /** basepr_[i][j] = base to which the jth base is paired in the ith structure */
  int **basepr_;
  int ndbl_;
  int dbl_[maxforce];
  /** The energy (the indices refer to the individual structures held in this object). 
   * 10 x the Gibb's free energy of the ith structure, this
   * is done so that integer math can be performed
   *  and the nearest tenth of a kcal/mol can be
   * followed*/
  int energy_[s_maxstructures+1];
  int inter_[3];
  int nnopair_;
  int nopair_[maxforce];
  int ngu_;
  int gu_[maxgu];
  /** Labels for the structures in this CT file.  a string of information for each of the structures*/
  std::string ctlabel_[s_maxstructures+1]; //[s_ctheaderlength];
  /** Nucleotides is a character array to store the sequence information -- 
   * this will allow the program to keep T and U from getting confused*/
  char *nucs_;
  bool intermolecular_;
  /** Has the class function allocate() beein called yet? */
  bool allocated_;
  bool templated_;
  bool **tem_;
  //int **fce;//[maxbases+1][2*maxbases]
 

 public:
  int createFromCTFile(const char * path);
  RNAStructure(const std::string &path);
  ~RNAStructure();
  int get_number_of_bases() const;
  int get_number_of_structures() const;
  std::string get_ith_label(int i) const;
  std::string get_dot_parens_structure(int i) const;
  bool intermolecular() const;
  int ** basepr() { return basepr_; }
  inline int numseq(int i) { return numseq_[i]; }
 private:
  void allocate(int size = s_maxbases);
  void allocatetem();
  void tonum(char *base, int count);
  
};


/**
 *  contains a stack of data, used by
 * functions that analyze a structure piecewise
 */
class Structstack {
  /** The actual stack. Indices=? */
  int stk_[51][4];
  /** stackpointer */
  int sp_;
 public:
  Structstack();
  void push(int a, int b, int c, int d);
  void pull(int *i, int *j, int *open, int *null, int *stz);
};




#endif
