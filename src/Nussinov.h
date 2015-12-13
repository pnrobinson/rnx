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
 * \author (last to touch it) $Author: bv $
 *
 * \version $Revision: 0.0.1 $
 *
 * \date $Date: 2005/04/14 14:16:20 $
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


  void traceback(int i, int j);
  
 public:
  Nussinov(const char * rna);
  ~Nussinov();
  /** Get length of the RNA sequence being analysed. */
  unsigned int get_len() const;
  /** Return a parenthesis-dot representation of the RNA secondary structure. */
  char * fold_rna();

  



};


#endif
