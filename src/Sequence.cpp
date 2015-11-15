

#include "Sequence.h"


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
