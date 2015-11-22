#ifndef SEQUENCE_H
#define SEQUENCE_H

/** 
 * Sequence.h
 * A set of classes designed to input and store DNA and RNA sequence data and
 * to provide a well defined interface to the RNA folding sections of the code.
 * @author Peter Robinson
 * @version 0.0.2 (22 November 2015)
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
  /** The NCBI gi number of the sequence */
  unsigned int gi_;
  /** The accession number plus version suffix */
  std::string accession_;
  /** The name of the locus */
  std::string locus_;
  /** Start position (one based) of coding sequence. */
  unsigned int CDS_startpos_;
  /** End position (one based) of coding sequence. */
  unsigned int CDS_endpos_;
 public:
  Record(std::string h);
  void appendSequenceLine(std::string line);
  unsigned int get_size();
  std::string substr(unsigned int start_pos, unsigned int length);
  std::string get_rna() const;
  int get_gi() const;
  std::string get_accession_number() const;
  std::string get_locus() const;
  unsigned int get_CDS_startpos() const;
  unsigned int get_CDS_endpos() const;
  unsigned int get_CDS_length() const;
  std::string get_5utr() const;
  //
  void set_locus(std::string loc);
  void set_accession(std::string acc);
  void set_gi(int);
  void set_CDS_startpos(int startpos);
  void set_CDS_endpos(int endpos);
   // some implementation details
  void appendSequenceFromGeneBankLines(std::vector<std::string> seqlines);
};

bool parseFASTA(std::string path, std::vector<Record> & records);
bool parseGenBank(std::string path, std::vector<Record> & records);
 

#endif
/* eof */





