#ifndef NUSSINOV_H
#define NUSSINOV_H

/** 
 * Nussinov.h
 * A class to perform Nussinov-based alignments of RNA secondary structure.
 * Still prototyping.
 * @author Peter Robinson
 * @version 0.0.1 (7 December 2015)
 */




/**
 * \class Nussinov
 *
 * \ingroup Folding
 * (Note, this needs exactly one \defgroup somewhere)
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
  
 public:
  Nussinov(const char * rna);
  ~Nussinov();
  /** Get length of the RNA sequence being analysed. */
  unsigned int get_len() const;
  /** Return a parenthesis-dot representation of the RNA secondary structure. */
  char * fold_rna() const;

  



};


#endif
