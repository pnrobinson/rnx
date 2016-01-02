



#include "RNAStructure.h"

#include <iostream>
#include <fstream>

#include <cstring>
#include <sstream>
#include <stdlib.h>     /* atof */
#include <math.h>       /* floor */
#include <cerrno>


RNAStructure::RNAStructure(const std::string &path) {
  int i;
  allocate();
  for (i=1;i<=s_maxstructures;i++) {
    energy_[i]=0;
    allocated_ = false;
  }
  nnopair_=0;
  npair_=0;
  ndbl_=0;
  intermolecular_ = false;
  ngu_ = 0;
  templated_ = false;
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
* Note that this class was adapted from a struct that was intended to deal with CT files. It is
* thus not primarily adapted to dot-parens or other structure representations. We will make a
* function to output a dot-paren representation from the CT representation.
 * @param path Path to the ct file.
 * @param rnalist Vector that will hold all of the structures in the current CT file we are parsing.
 */
int RNAStructure::createFromCTFile(const char * path) {
  int count, i, j;
  int linelength = 20;
  char base[2]; // base[0] will be one or ACGU, and base[1]='\0'
  char line[linelength], temp[50];
  std::string header;
  std::string lin;
  std::ifstream in;
  numofbases_=-1; /* flag that numbases is not initialised */
  in.open(path);
  if (!in.is_open())
    perror("error while opening CT file");
  /*
  in >> count; // this is the number of residues in the sequence. 
  j = 0;

  if (count == -100) { //this is a CCT formatted file:
    std::cerr << "[ERROR] Attempting to read a CCT-formated file \"" << path << "\", which is not supported."
	      << std::endl;
    exit(1);
  }
  // Reset the file
  in.close();
  in.open(path);
  */
  for (numofstructures_ = 1; numofstructures_<=s_maxstructures;numofstructures_++) {
    if (numofstructures_==s_maxstructures) {
      std::cerr << "[ERROR] Number of structures in CT file \""<< path << "\" exceeds maximum allowed.\n";
      exit(1);
    }
    // The first line should also be the header. Format: 
    // 99	dG = -31.70 [Initially -31.70] AAEU02002163 1/2268-2170
    // alternative format: 
    // 76   ENERGY = 0.1  RA7680
    if (getline(in,header)) {
      /** Some CT files have an empty last line. This is probably not standard-conform, 
       * but if we see this we assume that the data is over and we stop trying to input
       * more data.
       */
      if (header.empty() || header.length()==0) {
	std::cerr << "[WARNING] Empty header line (or last empty line) in CT file: \""
		  << header << "\" in file: "
		  << path << "\"\n";
	goto bailout;
      }
      const std::string delims(" \t"); /* space or tab */
      size_t first_space = header.find_first_of(delims);
      std::string num = header.substr(0,first_space );
      int n = atoi(num.c_str());
      if (numofbases_<0)
	numofbases_=n;
      else if (numofbases_ != n) {
	std::cerr << "[ERROR] Malformed CT file: Divergent number of bases (n=" << n
		  <<") but we previously found numofbases_="<<numofbases_ << std::endl;
	exit(1); /* This is a major error that means the input file is corrupt. Do not try to recover. */
      }
      size_t rest = header.find_first_not_of(delims,first_space+1);
      size_t last = header.find_last_not_of(delims);
      header = header.substr(rest, (last-rest+1));
    } else {
      if (errno != 0) { /* if errno==0, there was no error, and we are simply finished with the file. */
	perror("Could not process header line of CT file");
      }
      break;
    }
    /* when we get here, we are finished with the header */
    if (in.eof()) {
      numofstructures_--;
      return 1;
    }
    ctlabel_[numofstructures_] = header;
    // The following code gets the structural information.
    for (count=1; count<= numofbases_; count++)	{
      in >> temp; //ignore base number in ctfile
      in >> base; //read the base
      strcpy(base+1, "\0");
      nucs_[count]=base[0];
      tonum(base, count); //convert base to numeric
      if (numseq_[count]==5) {  /* 5 is a flag for intermolecular interactions, see erg3 */
	intermolecular_ = true;
	inter_[j] = count; /* the index of the jth base with intermolecular interactions. */
	j++;
      }
      in >> temp; //ignore numbering
      in >> temp; //ignore numbering
      in >> basepr_[numofstructures_][count]; //read base pairing info
      in >> hnumber_[count]; //read historical numbering
    }
    /* when we get here,  we have read in all of the lines for each base in the current structure. */
    char c = in.peek();  // peek character
    if ( c == EOF ) {
      goto bailout;
    } else if (c==10) { /* char value of 10 means new line. We need to remove this character
			   in order to get the next call to getline to work correctly! */
      in.get(c);
    }
  }
 bailout:
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

int RNAStructure::get_number_of_structures() const {
  return numofstructures_;
}

std::string RNAStructure::get_ith_label(int i) const {
  if (i>numofstructures_) {
    std::cerr << "[ERROR] Index of sought after label too high: " << i << std::endl;
    return NULL;
  }
  return ctlabel_[i];
}

bool RNAStructure::intermolecular() const {
  return intermolecular_;
}

/**
 * Note that we use the information in the basepr_ array to create
 * the dot-paren representation. basepr_[i][j] = base to which the jth base is paired in the ith structure.
 * Note that basepr_ uses one-based numbering. If basepr_[i][j]>j, then we are openng a base pair to
 * a base that comes later in the sequence, thus, we write '('. If basepr_[i][j]==0, there is no base
 * pair, and we write nothing. If basepr_[i][j]<j, then we are closing a base pair that was opened before,
 * and we write ')'.
 * @param i The index of the structure to be output (note: structures begin at i=1, not i=0).
 */
std::string RNAStructure::get_dot_parens_structure(int i) const {
  if (i>numofstructures_) {
    std::cerr << "[ERROR] Index of sought after label too high: " << i << std::endl;
    return NULL;
  }
  std::string s(numofbases_, '.');
  for (int j=1;j<= numofbases_; ++j) {
    if (basepr_[i][j] != 0 &&  basepr_[i][j] > j)
      s[j-1] = '(';
    else if  (basepr_[i][j] != 0 &&  basepr_[i][j] < j)
      s[j-1]=')';
  }
  return s;
}




/*********  Structstack ******************/

Structstack::Structstack() {
  sp_ = 0;  //set stack counter
}

/**
 * The parameters a,b,c,d corresponds to the four columns of the array stk_
 * @param a index i (see pull)
 * @param b index j (see pull)
 * @param c open (see pull)
 * @param d null (see pull)
 */
void Structstack::push(int a, int b, int c, int d) {
  sp_++;
  stk_[sp_][0]= a;
  stk_[sp_][1]= b;
  stk_[sp_][2]= c;
  stk_[sp_][3]= d;
}




void Structstack::pull(int *i, int *j, int *open, int *null, int *stz)
{
  if (sp_==0) {
    *stz = 1;
    return;
  } else {
    *stz = 0;
    *i = stk_[sp_][0];
    *j = stk_[sp_][1];
    *open= stk_[sp_][2];
    *null= stk_[sp_][3];
    sp_--;
  }
}

