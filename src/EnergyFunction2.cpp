



#include "EnergyFunction2.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>


/**
 * All the information read from the
 * thermodynamic data files.
 */



Datatable::Datatable(const char* directory_path) {
  size_t len = strlen(directory_path);
  dirpath_ = new char[len+1];
  strcpy(dirpath_,directory_path);
  input_data();
}

Datatable::~Datatable() {
  delete [] dirpath_;
}

const char* Datatable::get_data_dir() const {
  return dirpath_;
}


bool file_exists(const char * path) {
  FILE *check;
   //check that all the files exist with a C i/o function
  if ((check = fopen(path, "r")) == NULL) {
    std::cerr << path << " missing\n";
    return false;
  } else
    return true;
}


void Datatable::input_data() {
  std::string files[] = {"loop.dat",};
  int n_elem = sizeof(files)/sizeof(files[0]);
  for (unsigned int i=0;i<n_elem;++i) {
    std::stringstream ss;
    ss << dirpath_ << "/" << files[i];
    std::string s = ss.str();
    if (! file_exists(s.c_str())) {
      std::cerr << "Error, could not find file " << s << std::endl;
    }
  }
  
  /* std::ifstream infile("thefile.txt");
  std::string line;
  while (std::getline(infile, line))*/

}
