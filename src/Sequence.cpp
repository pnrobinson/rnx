

#include "Sequence.h"
#include <algorithm>
#include <cstdlib>
#include <sstream>


void Record::appendSequenceLine(std::string line){
  this->sequence_ += line;
}

unsigned int Record::get_size() {
  return this->sequence_.length();
}

std::string Record::substr(unsigned int start_pos, unsigned int length){
  return this->sequence_.substr(start_pos,length);
}

std::string Record::get_rna() const {
  std::string rna = this->sequence_;
  for (unsigned i=0;i<rna.size();++i) {
    if (rna[i]=='T')
      rna[i]='U';
  }
  return rna;
}

std::string Record::get_accession_number() const {
  if (!accession_.empty()) return this->accession_;
  else return "?";
}

int Record::get_gi() const {
  return gi_;
}

void Record::set_gi(int gi) {
  gi_=gi;
}


/**
 * Remove leading and trailing white space
 */
std::string trim(std::string& str)
{
    size_t first = str.find_first_not_of(' ');
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last-first+1));
}

/**
 * split a string on delimiters and return a vector.
 * @param str String to be split
 * @param delimiter character to split on
 */
std::vector<std::string> split(std::string str, char delimiter) {
  std::vector<std::string> internal;
  std::stringstream ss(str); // Turn the string into a stream.
  std::string tok;
  while(getline(ss, tok, delimiter)) {
    internal.push_back(trim(tok));
  }
  return internal;
}

void Record::set_locus(std::string loc){
  this->locus_ = loc;
}

std::string Record::get_locus() const{
  return this->locus_;
}
void Record::set_accession(std::string acc){
  this->accession_ = acc;
}


/**
 * Note that this constructor assumes we
 * are gettting an NCBI FASTA file. In that case, 
 * the header line often has a certain structure,
 * e.g., gi|21434723 for the NCBI "gi" number and
 * gb|accession|locus for GenBank sequences (e.g., gb|M73307|AGMA13GT).
 * @param h The header, e.g, >gi|21434723|gb|M73307|AGMA13GT| name...
 */
Record::Record(std::string h) : header_(h), sequence_() {
   /* try to parse the fields of the header. 
   */
  if (header_.length()>0 && header_.at(0)=='>') /* remove '>' if necessary */
    header_.replace(header_.begin(),header_.end(),header_.substr(1));
  
  char delimiter = '|';
  std::vector<std::string> vec = split(header_,delimiter);
  int len = vec.size();
  for (unsigned i=0;i<len;++i) {
    if (vec[i].compare("gi")==0){
      if (++i < len) {
	int gi = std::atoi (vec[i].c_str());
	gi_= gi;
      }
    } else if (vec[i].compare("gb")==0 || vec[i].compare("ref")==0) {
      if (++i < len) {
	accession_ = vec[i];
      }
      if (++i < len)
	locus_ = vec[i];
    }
  }
}

  /**
   * Parse a FASTA file and create one or more Record objects.
   */
bool parseFASTA(std::string path, std::vector<Record> & records){
  std::ifstream fin(path.c_str());
  if(!fin) {
    std::cerr << "Couldn't open the input FASTA file: \""
	      << path << "\"" << std::endl;
    return false;
  }

  std::string line;
  size_t counter = -1;
  // Priming read. May want to check that this line is a "header" before the loop.
  getline(fin, line);
  while(fin) {
    if(line[0] == '>')
    {
      records.push_back(Record(line));
      counter++;
    }
    else
    {
      records[counter].appendSequenceLine(line);
    }
    getline(fin, line);
  }
  fin.close();
 
  return true;
}

/**
 * This function removes leading space and digits as well
 * as any space, and returns the resulting string.
 * It is intended to deal with GenBank lines such as
 * "   601 taataaaaaa catttatttt cattgc"
 */
std::string remove_if_space_or_digit(std::string ori) {
  unsigned i = 0;
  std::string str = ori;
  while (std::isspace(str[i]))
    i++;
  if (std::isdigit(str[i]))
    while (std::isdigit(str[i]))
      i++;
  std::string dest = str.substr(i);
  std::string::iterator end_pos = std::remove(dest.begin(), dest.end(), ' ');
  dest.erase(end_pos, str.end());
  return dest;
}


/**
 * Adds nucleotide sequence data to the Record object from a list of 
 * lines from the Genbank file (between ORIGIN and //) that also
 * contain numbers and spaces. This function removes the number and
 * the spaces and converts all nucleotides to upper case. Note that
 * the function currently only accepts AaCcGgTt nucleotides.
 * @param seqlines a vector of GenBank formated lines representing the sequence.
 */
void Record::appendSequenceFromGeneBankLines(std::vector<std::string> seqlines){
  if (seqlines.size()==0)
    return;
 
  unsigned int sz = 60*seqlines.size(); // there are maximally 60 nucleotides per line
  this->sequence_.reserve(sz);
  unsigned int i=0;
  std::vector<std::string>::iterator it=seqlines.begin();
  for (; it != seqlines.end(); ++it) {
    for (unsigned int j=0;j<it->size();++j) {
      switch(it->at(j)) {
      case 'a':
      case 'A':
	sequence_.push_back('A');
	break;
      case 'c':
      case 'C':
	sequence_.push_back('C');
	break;
      case 'g':
      case 'G':
	sequence_.push_back('G');
	break;
      case 't':
      case 'T':
	sequence_.push_back('T');
	break;
      }
    }
  }
  //std::cout << "seq=" << sequence_ << std::endl << " lentgh=" << sequence_.length() << std::endl;
}

/**
 * Parse a GenBank formated file with one or more entries.
 * @param path Path to the GenBank file
 * @param records Reference to a vector that will be filled with one Record for each GenBank entry.
 */
bool parseGenBank(std::string path, std::vector<Record> & records) {
  std::ifstream fin(path.c_str());
  if(!fin) {
    std::cerr << "Couldn't open the input GenBank file: \""
	      << path << "\"" << std::endl;
    return false;
  }

  std::string line;
  size_t counter = -1;
  // Priming read. May want to check that this line is a "header" before the loop.
  getline(fin, line);
  
  while(fin) {
    if(line.find("LOCUS") ==0) /* line starts with LOCUS */
    {
      records.push_back(Record(line));
      counter++;
      //e.g., DEFINITION Homo sapiens hemoglobin, beta (HBB), mRNA.
    } else if (line.find("DEFINITION")==0) {
      unsigned i=10;
      while (std::isspace(line.at(i)))
	i++;
      // Assume if we get here that there must be a Record object in the vector
      records.back().set_locus(line.substr(i));
    } else if (line.find("ACCESSION")==0) {
      unsigned i=10;
      while (std::isspace(line.at(i)))
	i++;
      records.back().set_accession(line.substr(i));
    } else if (line.find("VERSION")==0) {
       std::size_t found = line.find("GI:");
       if (found!=std::string::npos && found < 3+line.length()) {
	 std::string gi = line.substr(found+3);
	 int g = atoi(gi.c_str());
	 records.back().set_gi(g);
       }
     
    } else if  (line.find("ORIGIN")==0) {
      /* i.e., we are at beginning of sequence block
ORIGIN      
        1 acatttgctt ctgacacaac tgtgttcact agcaacctca aacagacacc atggtgcatc
      (...)
      601 taataaaaaa catttatttt cattgc
//
*/
      std::vector<std::string> seqlines;
      getline(fin, line); // advance to first sequence line
      while (fin) {
	if (line.find("//") == 0) {
	  break;
	} else {
	  seqlines.push_back(line);
	}
	getline(fin, line);
      }
      records.back().appendSequenceFromGeneBankLines(seqlines);
    }
      /* when we get here we should have just read the // line at the end of the sequence */
   /* else if ORIGIN */
    getline(fin, line);
  }
  fin.close();
  
  return true;
}




/* eof */
