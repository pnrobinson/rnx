#ifndef NUSSINOV_H
#define NUSSINOV_H




/**
 * \defgroup Folding RNA Folding Algorithms
 */


/**
 * \class Nussinov
 *
 * \ingroup Folding
 * 
 *
 * \brief Class for folding an RNA sequence according to Nussinov
 *
 * Right now, FASTA
 *
 * \note Still prototyping.
 *
 * \author Peter Robinson
 *
 * \version  0.0.2
 *
 * \date 19 December 2015
 *
 * Contact: peter.robinson@charite.de
 *
 * Created on: 9 December 2015
 *
 *
 */

class Nussinov {
  /** A copy of the RNA sequence we are to investigate (note, maybe we do not need to make copy). */
  char * rna_;
  /** This is the minimum distance between two nucleotides that undergo base pairing. The
      default is set to 3.*/
  unsigned int h_;
  /** The Nussinov matrix containing the number of bonds for the subsequence i,j) */
  int ** mat_;
  /** The parenthesis - dot representation of the secondary structure of the current RNA */
  char * structure_;


  /** Find an optimal secondary structure (this method is called by fold_rna). */
  void traceback();
  /** Can be used to print out the DP matrix for debugging purposes. */
  void debugPrintMatrix();

  void init(const char * rna);
  
 public:
  Nussinov(const char * rna);
  Nussinov(const char * rna, unsigned int h);
  ~Nussinov();
  unsigned int get_len() const;
  const char * fold_rna();
  
};


#endif
