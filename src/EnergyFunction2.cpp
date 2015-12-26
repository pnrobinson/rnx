



#include "EnergyFunction2.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <stdlib.h>     /* atof */
#include <math.h>       /* floor */

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


/**
 * This function performs a quick test of file existence
 * by opening the file pointed to be path (or attempting to).
 * If the fopen returns NULL, the file is not there, return false.
 */
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
      //exit(1);
    }
  }
  /* If we get here, all of the files we need can be found. */
  std::stringstream ss;
  ss << dirpath_ << "/loop.dat";
  std::string s = ss.str();
  input_loop_dat(s);
  
  
   
}

/**
 * Input the file "loop.dat" (referred to as loop2/ lo1.open(loop2);  in e2f2).
 * loop.dat has information about  internal loops, hairpin loops, and bulge loops 
 * Note that we first skip down past the line that starts with "-----" ...
 * @param path Path to loop.dat
 /*	read info from loop for
  for (i=1; i <= 30; i++) {
    lo1 >> lineoftext; //get past the size column in table
    lo1 >> lineoftext;
    if (strcmp(lineoftext, ".")){
      data->inter[i] = (int) floor (100.0*(atof(lineoftext))+.5);
    }
    else data->inter[i] = infinity;

    //cout <<"inter = "<<data->inter[i]<<"\n";
    lo1 >> lineoftext;
    if (strcmp(lineoftext, "."))
      data->bulge[i] = (int) floor(100.0*(atof(lineoftext))+.5);
    else data->bulge[i] = infinity;

    //cout <<"bulge = "<<data->bulge[i]<<"\n";
    lo1 >> lineoftext;
    if (strcmp(lineoftext, ".")){
      data->hairpin[i] = (int) floor(100.0*(atof(lineoftext))+.5);
    }
    else data->hairpin[i] = infinity;

    //cout <<"hair = "<<data->hairpin[i]<<"\n";
  }
 */
void Datatable::input_loop_dat(std::string &path){
  std::ifstream infile(path.c_str());
  std::string token;
  int c=0;
  while ( infile >> token)  {
    if (token.substr(0,5) == "-----")
      break;
  }
  for (unsigned int i=1; i <= 30; i++) {
    infile >> token; //get past the size column in table
    infile >> token;
    if (token == "."){
      inter_[i] = s_infinity;
    } else {
      inter_[i] = static_cast<int> (floor (100.0*(atof(token.c_str()))+.5)) ;
    }
    //std::cout << i << ")" << inter[i] << std::endl;
    infile >> token;
    if (token ==  ".") {
      bulge_[i] = s_infinity;
    } else {
      bulge_[i] = (int) floor(100.0*(atof(token.c_str()))+.5);
    }
    //std::cout <<"bulge = "<<data->bulge[i]<<"\n";
    infile >> token;
    if (token ==  ".") {
      hairpin_[i] = s_infinity;
    } else {
      hairpin_[i] = (int) floor(100.0*(atof(token.c_str()))+.5);
    }
  }
}




int Datatable::get_destabilizing_energy_internal_loop(unsigned int loop_size) {
  return inter_[loop_size];
}

int Datatable::get_destabilizing_energy_bulge_loop(unsigned int loop_size) {
  return bulge_[loop_size];
}

int Datatable::get_destabilizing_energy_hairpin_loop(unsigned int loop_size) {
  return hairpin_[loop_size];
}
