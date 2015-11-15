#ifndef SEQUENCE_H
#define SEQUENCE_H

/** ===========================================================================
 * Sequence.h
 * A set of classes designed to input and store DNA and RNA sequence data and
 * to provide a well defined interface to the RNA folding sections of the code.
 * @author Peter Robinson
 * @version 0.0.1 (15 November 2015)
 * ============================================================================
 */




#include <iostream>
#include <string>
#include <fstream>
#include <vector>


/**
 * \class Record
 *
 * \ingroup PackageName
 * (Note, this needs exactly one \defgroup somewhere)
 *
 * \brief Class for representing a single sequence from a file
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
 * Created on: Wed Apr 13 18:39:37 2005
 *
 * $Id: doxygen-howto.html,v 1.5 2005/04/14 14:16:20 bv Exp $
 *
 */
class Record {
 private:
  /** The name */
  std::string header_;
  /** The DNA or RNA sequence */
  std::string sequence_;
 public:
 Record(std::string h) : header_(h), sequence_() {}
  void appendSequenceLine(std::string line);
  unsigned int get_size();
  std::string substr(unsigned int start_pos, unsigned int length);
  std::string get_rna() const;

};

bool parseFASTA(std::string path, std::vector<Record> & records);
 

#endif
/* eof */





