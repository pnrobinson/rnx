



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

/**
 * @return the directory in which the thermodynamic data files are stored
 */
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
  std::string files[] = {"loop.dat","stack.dat", "tstackh.dat",
			 "tstacki.dat", "tloop.dat", "miscloop.dat"};
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
  /* 1) loop */
  std::stringstream ss;
  ss << dirpath_ << "/loop.dat";
  std::string s = ss.str();
  input_loop_dat(s);
  /* 2) stack */
  ss.str(std::string()); /* reset */
  ss << dirpath_ <<  "/stack.dat";
  s = ss.str();
  input_stack_dat(s);
  /* 3) tstackh */
  ss.str(std::string()); /* reset */
  ss << dirpath_ <<  "/tstackh.dat";
  s = ss.str();
  input_tstackh_dat(s);
  /* 4) tstacki */
  ss.str(std::string()); /* reset */
  ss << dirpath_ <<  "/tstacki.dat";
  s = ss.str();
  input_tstacki_dat(s);
  /* 5) tloop */
  ss.str(std::string()); /* reset */
  ss << dirpath_ <<  "/tloop.dat";
  s = ss.str();
  input_tloop_dat(s);
   /* 6) miscloop */
  ss.str(std::string()); /* reset */
  ss << dirpath_ <<  "/miscloop.dat";
  s = ss.str();
  input_miscloop_dat(s);
   
}

/**
 * Input the file "loop.dat" (referred to as loop2/lo1 in e2f2).
 * loop.dat has information about  internal loops, hairpin loops, and bulge loops 
 * Note that we first skip down past the line that starts with "-----" ...
 * then there are four columns with index/internal/bulge/hairpin with 30 lines.
 * @param path Path to loop.dat
 */
void Datatable::input_loop_dat(std::string &path){
  std::ifstream infile(path.c_str());
  if (! infile.is_open() ) {
    std::cerr << "[ERROR] could not initialize " << path << " for I/O" << std::endl;
  }
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


/**
 * Read info from stack 
 * add to the stack table the case where X (represented as 0) is looked up:
 * \verbatim
   Data Arangement: 
 
                 Y           
         ------------------  
     (X)  A    C    G    U   
         ------------------  
             5' ==> 3'       
                AX           
                AY 
             3' <== 5'       
     (A)   .     .     .     .   
     (C)   .     .     .     .   
     (G)   .     .     .     .   
     (U) -0.7  -0.1  -0.7  -0.1  
  \endverbatim
  * The total free energy of the entire Stacking Region is given by 
  * the addition of each pair of  adjacent base pairs, including 
  * energy contributions for both base pair stacking and hydrogen bonding.
  * The Stacking Energies table being parsed here  uses the following
  * arrangement for a stack:
  \verbatim 5’–WX–3’
  3’–ZY–5’
  \endverbatim
  * The corresponding energy would appear in the Wth row and the Zth
  * column of 4 by 4  tables, and in the Xthrow and the Yth column of that table.
  * The file <b>stack.dat</b> has four rows and four columns. If we index this
  * arrangement with w,x, and if we index the rows and columns of the individual
  * tables with y, z, then stack_[w][x][y][z] will give the correct values (see above). Note
  * that 1-based numbering is being used.
  *
  */
void Datatable::input_stack_dat(std::string &path) {
  std::ifstream infile(path.c_str());
  if (! infile.is_open() ) {
    std::cerr << "[ERROR] could not initialize " << path << " for I/O" << std::endl;
  }
  std::string token;
  int count, i, j, k , l;

  for (count=1; count<=42; count++)
    infile >> token; //get past text in file
  for (i=0; i<=5; i++) {
    if ((i!=0)&&(i!=5))
      for (count=1; count<=60; count++)
	infile >> token;
    for (k=0; k<=5; k++) {
      for (j=0; j<=5; j++) {
	for (l=0; l<=5; l++) {
	  if ((i==0)||(j==0)||(k==0)||(l==0)) {
	    stack_[i][j][k][l]=0;
	  }
	  else if ((i==5)||(j==5)||(k==5)||(l==5)) {
	    stack_[i][j][k][l] = s_infinity;
	  } else {
	    infile >> token;
	    if (token != "."){
	      stack_[i][j][k][l] =(int) floor(100.0*(atof(token.c_str()))+.5);
	    } else {
	      stack_[i][j][k][l] = s_infinity;
	    }
	  }
	}
      }
    }
  }
}


/**
 * Read info from tstackh 
 * add to the tstackh table the case where X (represented as 0) is looked up.
 * Function and datastructures analogous to input_stack_dat.
*/
void Datatable::input_tstackh_dat(std::string &path) {
  std::ifstream infile(path.c_str());
  if (! infile.is_open() ) {
    std::cerr << "[ERROR] could not initialize " << path << " for I/O" << std::endl;
  }
  std::string token;
  int count, i, j, k , l;

  for (count=1; count<=46; count++)
    infile >> token; //get past text in file
  for (i=0; i<=5; i++) {
    if ((i!=0)&&(i!=5))
      for (count=1; count<=60; count++)
	infile >> token;
    for (k=0; k<=5; k++) {
      for (j=0; j<=5; j++) {
	for (l=0; l<=5; l++) {
	  if ((i==0)||(j==0)||(k==0)||(l==0)) {
	    tstkh_[i][j][k][l]=0;
	  } else if ((i==5)||(j==5)) {
	    tstkh_[i][j][k][l] = s_infinity;
	  } else if ((i!=5)&&(j!=5)&&((k==5)||(l==5))) {
	    tstkh_[i][j][k][l] = 0;
	  } else {
	    infile >> token;
	    if (token !=  ".") {
	      tstkh_[i][j][k][l] = static_cast<int> (floor(100.0*(atof(token.c_str()))+.5));
	    } else {
	      tstkh_[i][j][k][l] = s_infinity;
	    }
	  }
	}
      }
    }
  }
}


/**
 * Read info from tstacki 
 * add to the tstacki table the case where X (represented as 0) is looked up.
 * Function and datastructures analogous to input_stack_dat.
*/
void Datatable::input_tstacki_dat(std::string &path) {
  std::ifstream infile(path.c_str());
  if (! infile.is_open() ) {
    std::cerr << "[ERROR] could not initialize " << path << " for I/O" << std::endl;
  }
  std::string token;
  int count, i, j, k , l;

  for (count=1; count<=46; count++)
    infile >> token; //get past text in file
  for (i=0; i<=5; i++) {
    if ((i!=0)&&(i!=5))
      for (count=1; count<=60; count++)
	infile >> token;
    for (k=0; k<=5; k++) {
      for (j=0; j<=5; j++) {
	for (l=0; l<=5; l++) {
	  if ((i==0)||(j==0)||(k==0)||(l==0)) {
	    tstki_[i][j][k][l]=0;
	  } else if ((i==5)||(j==5)) {
	    tstki_[i][j][k][l] = s_infinity;
	  } else if ((i!=5)&&(j!=5)&&((k==5)||(l==5))) {
	    tstki_[i][j][k][l] = 0;
	  } else {
	    infile>> token;
	    if (token != "."){
	      tstki_[i][j][k][l]=(int)floor(100.0*(atof(token.c_str()))+.5);
	    } else {
	      tstki_[i][j][k][l] = s_infinity;
	    }
	  }
	}
      }
    }
  }
}



void Datatable::input_miscloop_dat(std::string &path) {

}


/**
 * Used to create a numerical code for sequences in input_tloop_dat
 */
int tonumi(const char *base)	{
  int	a;
  if (!strcmp(base, "A")||!strcmp(base, "B")) (a = 1);
  else if(!strcmp(base, "C")||!strcmp(base, "Z")) (a = 2);
  else if(!strcmp(base, "G")||!strcmp(base, "H")) (a = 3);
  else if(!strcmp(base, "U")||!strcmp(base, "V")) (a = 4);
  else if(!strcmp(base, "T")||!strcmp(base, "W")) (a = 4);
  else (a=0);  //this is for others, like X
  return a;
}

/**
 * Input file "tloop.dat". The tetraloop bonus table. For hairpin loops with four unpaired nucleotides, this table is consulted. 
 * If the sequence (starting with the 5' last paired nucleotide and finishing with its 3' paired nucleotide) appears in the table, 
 * the bonus is applied to the hairpin's stability.
*/
void Datatable::input_tloop_dat(std::string &path) {
  std::ifstream infile(path.c_str());
  char base[110];
  char temp[300];
  if (! infile.is_open() ) {
    std::cerr << "[ERROR] could not initialize " << path << " for I/O" << std::endl;
  }
  std::string token;
  int count;
 
   /*	Read info from tloops */
  for (count=1; count<=3; count++)
    infile >> token; //get past text in file
  numoftloops_=0;
  infile >> token;


  for (count=1; count<=s_maxtloop&&!infile.eof(); count++){
    //cout << lineoftext;
    numoftloops_++;
    strcpy(base, token.c_str());
    strcpy(base+1, "\0");
    tloop_[numoftloops_][0] = tonumi(base);
    //std::cout << base << " ";
    //std::cout << tloop_[numoftloops_][0] << "\n";

    //strcpy(base, token[1]);
    base[0] = token[1];
    strcpy(base+1, "\0");
    tloop_[numoftloops_][0] = tloop_[numoftloops_][0]+  5*tonumi(base);
    //std::cout << base << " ";
    //std::cout << tloop_[numoftloops_][0] << "\n";

    //strcpy(base, token+2);
    base[0] = token[2];
    strcpy(base+1, "\0");
    tloop_[numoftloops_][0] = tloop_[numoftloops_][0]+  25*tonumi(base);
    //cout << base << "\n";
    //cout << data->tloop[data->numoftloops][0] << "\n";
    //strcpy(base, token+3);
    base[0] = token[3];
    strcpy(base+1, "\0");
    tloop_[numoftloops_][0] = tloop_[numoftloops_][0]+  125*tonumi(base);
    //cout << base << "\n";
    //cout << data->tloop[data->numoftloops][0] << "\n";
    //strcpy(base, lineoftext+4);
    base[0] = token[4];
    strcpy(base+1, "\0");
    tloop_[numoftloops_][0] = tloop_[numoftloops_][0]+  625*tonumi(base);
    //cout << base << "\n";
    //cout << data->tloop[data->numoftloops][0] << "\n";
    //strcpy(base, lineoftext+5);
    base[0] = token[5];
    strcpy(base+1, "\0");
    tloop_[numoftloops_][0] = tloop_[numoftloops_][0]+ 3125*tonumi(base);
    //cout << base << "\n";
    //cout << data->tloop[data->numoftloops][0] << "\n";
    infile >> temp;
    tloop_[numoftloops_][1] = (int) floor (100.0*atof(temp)+0.5);
    //std::cout<< "ntloops[0]="<< numoftloops_ << ": "<<base <<  ":"<< tloop_[numoftloops_][0] << std::endl;
    //std::cout<< "ntloops[1]="<< numoftloops_ << ": "<<  tloop_[numoftloops_][1] << std::endl;
    infile >> token;
  }

}

int Datatable::get_tetraloop_energy(const char *seq) const {
  int y=0;
  size_t len = strlen(seq);
  if (len != 6)
    return 0;
  y = tonumi(seq);
  y = y +  5*tonumi(seq+1);
  y = y +  25*tonumi(seq+2);
  y = y +  125*tonumi(seq+3);
  y = y +  625*tonumi(seq+4);
  y = y +  3125*tonumi(seq+5);
  std::cout << "get tloop y=" << y << std::endl;
  for (int i=1;i<s_maxtloop;++i) {
    std::cout << "tloop[" << i << "]=" << tloop_[i][0] << std::endl;
    if (tloop_[i][0] == y)
      return tloop_[i][1];
  }
  return 0;
}



int Datatable::get_destabilizing_energy_internal_loop(unsigned int loop_size) const {
  return inter_[loop_size];
}

int Datatable::get_destabilizing_energy_bulge_loop(unsigned int loop_size) const {
  return bulge_[loop_size];
}

int Datatable::get_destabilizing_energy_hairpin_loop(unsigned int loop_size) const {
  return hairpin_[loop_size];
}

int Datatable::get_stack_energy(int w, int x, int y, int z) const {
  return stack_[w][x][y][z];
}

/**
 * @return terminal stacking energy of base-pair or mismatch 
 */
int Datatable::get_tstackh_energy(int w, int x, int y, int z) const {
  return tstkh_[w][x][y][z];
}

/**
 * @return terminal stacking energy of base-pair or mismatch 
 */
int Datatable::get_tstacki_energy(int w, int x, int y, int z) const {
  return tstki_[w][x][y][z];
}



