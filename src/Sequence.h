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

};


std::vector<Record> parseFASTA(std::string path);
 

#endif
/* eof */





