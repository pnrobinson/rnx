

#include "Sequence.h"
#include <sstream>
#include <cstdlib>

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
 * \TODO Make robust GenBank parser!
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
    } else {
      records[counter].appendSequenceLine(line);
    }
    getline(fin, line);
  }
  fin.close();
 
  return true;
}
