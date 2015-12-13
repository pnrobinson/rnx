#ifndef NUSSINOV_H
#define NUSSINOV_H

/** 
 * Nussinov.h
 * A class to perform Nussinov-based alignments of RNA secondary structure.
 * Still prototyping.
 * @author Peter Robinson
 * @version 0.0.2 (11 December 2015)
 */


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
 * \note Attempts at zen rarely work.
 *
 * \author Peter Robinson
 *
 * \version  0.0.2
 *
 * \date $Date: 2015/12/11 $
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
  /** The Nussinov matrix containing the number of bonds for the subsequence i,j) */
  int ** mat_;
  /** The parenthesis - dot representation of the secondary structure of the current RNA */
  char * structure_;


  void traceback2(int i, int j);
  void traceback();
  void debugPrintMatrix();
  
 public:
  /** @param rna a string (upper case) of RNA nucleotides */
  Nussinov(const char * rna);
  ~Nussinov();
  /** Get length of the RNA sequence being analysed. */
  unsigned int get_len() const;
  /** Return a parenthesis-dot representation of the RNA secondary structure. Note that clients should make a copy of this string if they need to alter it*/
  const char * fold_rna();

  



};


#endif
