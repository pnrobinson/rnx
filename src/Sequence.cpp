

#include "Sequence.h"


void Record::appendSequenceLine(std::string line){
  this->sequence_ += line;
}

unsigned int Record::get_size() {
  return this->sequence_.length();
}

std::vector<Record> parseFASTA(std::string path) {
  std::vector<Record> records;
  std::ifstream fin(path.c_str());
  if(!fin) {
    std::cerr << "Couldn't open the input FASTA file: \""
	      << path << "\"" << std::endl;
    return records;
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
  
  return records;
}
