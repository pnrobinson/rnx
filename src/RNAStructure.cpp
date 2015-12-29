



#include "RNAStructure.h"

#include <iostream>
#include <fstream>

#include <cstring>
#include <sstream>
#include <stdlib.h>     /* atof */
#include <math.h>       /* floor */



RNAStructure::RNAStructure(const std::string &path) {
  createFromCTFile(path.c_str());
}

/**
 * The connect format is column based. The first column specified the sequence index, starting at one. 
 * Columns 3, 4, and 6 redundantly give sequence indices (plus/minus one). The second column contains 
 * the base in one-letter notation. Column 4 specifies the pairing partner of this base if it involved 
 * in a base pair. If the base is unpaired, this column is zero.
 * The parser expects one header line containing the word "ENERGY", "Energy", or "dG". Arbitrary header 
 * lines will not be accepted. Files in connect format may contain multiple sequence/structure records. 
 * The specified pseudoknot removal method will be applied to all records in the file.
\verbatim
73 ENERGY =     -17.50    S.cerevisiae_tRNA-PHE
 1 G       0    2   72    1
 2 C       1    3   71    2
 3 G       2    4   70    3
 4 G       3    5   69    4
 5 A       4    6   68    5
 6 U       5    7   67    6
 7 U       6    8   66    7
 8 U       7    9    0    8

              .
              .
              .

66 A      65   67    7   66
67 A      66   68    6   67
68 U      67   69    5   68
69 U      68   70    4   69
70 C      69   71    3   70
71 G      70   72    2   71
72 C      71   73    1   72
73 A      72   74    0   73
\endverbatim
* Note that the example file  RA7680.ct (used for unit testing and saved in the test_data directory)
* represents a tRNA (Sprinzl M, Vassilenko KS. Compilation of tRNA sequences and sequences of 
* tRNA genes. Nucleic Acids Res. 2005;33:D139–140).
 * @param path Path to the ct file.
 * @param rnalist Vector that will hold all of the structures in the current CT file we are parsing.
 */
int RNAStructure::createFromCTFile(const char * path) {
  int count, i, j;
  int linelength = 20;
  //RNAStructure rnastruct;
  //int ctheaderlength = 125; //max length of string containing info on sequence
  char base[2], header[s_ctheaderlength], line[linelength], temp[50];
  std::ifstream in;
  in.open(path);
  in >> count; /* this is the number of residues in the sequence. */
  j = 0;

  if (count == -100) { //this is a CCT formatted file:
    std::cerr << "[ERROR] Attempting to read a CCT-formated file \"" << path << "\", which is not supported."
	      << std::endl;
    exit(1);
  }
  // allocate memory
  allocate(count);
  // Rest the file
  in.close();
  in.open(path);
  for (numofstructures_ = 1; numofstructures_<=s_maxstructures;
       numofstructures_++) {
    strcpy (header, "");
    if (numofstructures_==s_maxstructures) {
      std::cerr << "[ERROR] Number of structures in CT file \""<< path << "\" exceeds maximum allowed.\n";
      exit(1);
    }
    in >> numofbases_;
    strcpy(line, "");
    in.getline(header, s_ctheaderlength);

      //	 do {
      //	 	strcat(header, line);
      //		if (in.eof()) {
      //			ct->numofstructures--;
      //			return;
      //		}
      //		in >> line;
      //		strcat(header, " ");
      //	 }
      //	 while(strcmp(line, "1"));
      //	 in.putback(*line); */

    if(in.eof()) {
      numofstructures_--;
      return 1;
    }
    
    strcpy((ctlabel_[numofstructures_]), header);
    for (count=1; count<= numofbases_; count++)	{
      
      //if(in.eof()) {
      //	ct->numofstructures--;
      //	return;
      //}
      in >> temp; //ignore base number in ctfile
      in >> base; //read the base
      strcpy(base+1, "\0");
      nucs_[count]=base[0];
      tonum(base, count); //convert base to numeric
      if (numseq_[count]==5) {
	intermolecular_ = true;
	inter_[j] = count;
	j++;
      }
      in >>	temp; //ignore numbering
      in >> temp; //ignore numbering
      in >> basepr_[numofstructures_][count]; //read base pairing info
      in >> hnumber_[count]; //read historical numbering
    }
  }
  numofstructures_--;
  return 1;
}



void RNAStructure::tonum(char *base, int count)	{
  if (!strcmp(base, "A")) (numseq_[count] = 1);
  else if(!strcmp(base, "B")) {
    (numseq_[count] = 1);

  }
  else if(!strcmp(base, "a")) {
    numseq_[count]=1;
    nnopair_++;
    nopair_[nnopair_] = count;
  }
  else if(!strcmp(base, "C")) (numseq_[count] = 2);
  else if(!strcmp(base, "Z")) {
    (numseq_[count] = 2);

  }
  else if(!strcmp(base, "c")) {
    numseq_[count] = 2;
    nnopair_++;
    nopair_[nnopair_] = count;
  }
  else if(!strcmp(base, "G")) (numseq_[count] = 3);
  else if(!strcmp(base, "H")) {
    (numseq_[count] = 3);

  }
  else if(!strcmp(base, "g")) {
    numseq_[count] = 3;
    nnopair_++;
    nopair_[nnopair_] = count;
  }

  else if(!strcmp(base, "U")||!strcmp(base, "T")) (numseq_[count] = 4);
  else if(!strcmp(base, "V")||!strcmp(base, "W")) {
    (numseq_[count] = 4);

  }
  else if(!strcmp(base, "u")||!strcmp(base, "t")) {
    numseq_[count] = 4;
    nnopair_++;
    nopair_[nnopair_] = count;
  }

  else if(!strcmp(base, "I")) {
    numseq_[count] = 5;
    intermolecular_= true;
  }

  else (numseq_[count]=0);  //this is for others, like X
  return;
}


void RNAStructure::allocate(int size) {
  int i;
  //Size = size;//save the size of the array so that the destructor can
  //deallocate the space
  numseq_ = new int [2*size+1];
  hnumber_ = new int [size+1];
  nucs_ = new char [size+2];
  basepr_ = new int *[maxstructures+1];
  for (i=0;i<=s_maxstructures;i++) {
    basepr_[i] = new int [size+1];
  }
  allocated_ = true;
}

RNAStructure::~RNAStructure() {
  int i;
  if (allocated_) {
    delete[] numseq_;
    for (i=0;i<=s_maxstructures;i++) {
      delete[] basepr_[i];
    }
    delete[] basepr_;
    delete[] hnumber_;
    delete[] nucs_;
  }
  if (templated_) {
    for (i=0;i<=numofbases_;i++) {
      delete[] tem_[i];
    }
    delete[] tem_;
  }
}



int RNAStructure::get_number_of_bases() const {
  return numofbases_;
}
