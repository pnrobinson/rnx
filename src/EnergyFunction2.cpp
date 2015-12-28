



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
  if (directory_path[len-1] == '/') {
    len--; /* Only copy up to but not including the final '/' if there is one */
  }
  dirpath_ = new char[len+1];
  strncpy(dirpath_,directory_path,len);
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
			 "tstacki.dat", "tloop.dat", "miscloop.dat",
			 "dangle.dat","int22.dat","int21.dat",
			 "coaxial.dat","triloop.dat", "tstackcoax.dat",
			 "coaxstack.dat","tstack.dat","tstackm.dat",
			 "int11.dat",};
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
  /* 7) dangle */
  ss.str(std::string()); /* reset */
  ss << dirpath_ <<  "/dangle.dat";
  s = ss.str();
  input_dangle_dat(s);
  /* 8) int22 */
  ss.str(std::string()); /* reset */
  ss << dirpath_ <<  "/int22.dat";
  s = ss.str();
  input_int22_dat(s);
  /* 9) int21 */
  ss.str(std::string()); /* reset */
  ss << dirpath_ <<  "/int21.dat";
  s = ss.str();
  input_int21_dat(s);
  /* 10) coaxial */
  ss.str(std::string()); /* reset */
  ss << dirpath_ <<  "/coaxial.dat";
  s = ss.str();
  input_coaxial_dat(s);
  /* 10) tstackcoax */
  ss.str(std::string()); /* reset */
  ss << dirpath_ <<  "/tstackcoax.dat";
  s = ss.str();
  input_tstackcoax_dat(s);
  /* 11) coaxstack */
  ss.str(std::string()); /* reset */
  ss << dirpath_ <<  "/coaxstack.dat";
  s = ss.str();
  input_coaxstack_dat(s);
  /* 12) tstack */
  ss.str(std::string()); /* reset */
  ss << dirpath_ <<  "/tstack.dat";
  s = ss.str();
  input_tstack_dat(s);
  /* 12) tstackm */
  ss.str(std::string()); /* reset */
  ss << dirpath_ <<  "/tstackm.dat";
  s = ss.str();
  input_tstackm_dat(s);
  /* 11) int11 */
  ss.str(std::string()); /* reset */
  ss << dirpath_ <<  "/int11.dat";
  s = ss.str();
  input_int11_dat(s);
   
}

/**
 * Input the file "loop.dat" (referred to as loop2/lo1 in e2f2).
 * loop.dat has information about  internal loops, hairpin loops, and bulge loops 
 * Note that we first skip down past the line that starts with "-----" ...
 * then there are four columns with index/internal/bulge/hairpin with 30 lines.
 * @param path Path to loop.dat
 */
void Datatable::input_loop_dat(const std::string &path){
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
void Datatable::input_stack_dat(const std::string &path) {
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
void Datatable::input_tstackh_dat(const std::string &path) {
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
void Datatable::input_tstacki_dat(const std::string &path) {
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


/**
 * Input the file "miscloop.dat" (Miscellaneous free energy rules ).
 * the key sequence "-->" now indicates a record in that file
 */
void Datatable::input_miscloop_dat(const std::string &path) {
  std::ifstream infile(path.c_str());
  if (! infile.is_open() ) {
    std::cerr << "[ERROR] could not initialize " << path << " for I/O" << std::endl;
  }
  std::string token;
  std::string temp;
  int count, i, j, k , l;


  infile >> token;
  while(token !=  "-->")
    infile >> token;

  infile >> temp;
  prelog_ = atof(temp.c_str()) * 100.0;
  //data->prelog = (data->prelog)*100.0;

 
  infile >> token;
  while(token != "-->")
    infile >> token;
  infile >> temp;
  maxpen_ = static_cast<int>( atof(temp.c_str())*100.0 + .5);


  infile >> token;
  while(token != "-->")
    infile >> token;

  for (count=1; count<= 4; count ++){
    infile >> temp;
    poppen_[count] = static_cast<int> (atof(temp.c_str())*100.0 + .5);
  } 										//this reads float values, converts
  // 	them int and assigns them into
  //		array poppen

  
  infile >> token;
  while(token != "-->")
    infile >> token;
  // assign some variables that are"hard-wired" into code
  eparam_[1] = 0; 						
  eparam_[2] = 0; 					
  eparam_[3] = 0;
  eparam_[4] = 0;
  
  infile >> temp;
  eparam_[5] = static_cast<int> ( floor (atof(temp.c_str())*100.0+.5) );  //constant multi-loop penalty

  infile >> temp;
  eparam_[6] = static_cast<int> ( floor (atof(temp.c_str())*100.0+.5) ); //(int) floor (temp*100.0+.5);

  eparam_[7] = 30;
  eparam_[8] = 30;
  eparam_[9] = -500;

  infile >> temp;
  eparam_[10] =  static_cast<int> ( floor (atof(temp.c_str())*100.0+.5) ); //(int) floor (temp*100.0+.5);


  infile >> token;
  while(token != "-->")
    infile >> token;

  infile >> temp;

  if (infile.peek()==EOF) {
    //these are old energy rules -- treat the other constants properly
    std::cerr << "[WARNING]: Attempt to use old energy rules (miscloop.dat). "
	      << "Please obtain modern file that continues with \"efn2 multibranched loops\"" << std::endl;
    exit(1);
  }
 
  efn2a_ = static_cast<int> ( floor (atof(temp.c_str())*100.0+.5));  //constant multi-loop penalty for efn2

  infile >> temp;
  efn2b_=  static_cast<int> ( floor (atof(temp.c_str())*100.0+.5));  //(int) floor(temp*100.0+.5);

  infile >> temp;
  efn2c_=  static_cast<int> ( floor (atof(temp.c_str())*100.0+.5));  //(int) floor(temp*100.0+.5);

  //now read the terminal AU penalty:
  infile >> token;
  while(token != "-->")
    infile >> token;
  infile >> temp;

  auend_ = static_cast<int> ( floor (atof(temp.c_str())*100.0+.5)); //(int) floor (temp*100.0+.5);

  //now read the GGG hairpin bonus:
  infile>>token;
  while(token != "-->")
    infile>>token;
  infile >> temp;
  gubonus_ = static_cast<int> ( floor (atof(temp.c_str())*100.0+.5)); //(int) floor (temp*100.0+.5);
 
  //now read the poly c hairpin penalty slope:
  infile >> token;
  while(token != "-->")
    infile>>token;
  infile >> temp;
  cslope_ = static_cast<int> ( floor (atof(temp.c_str())*100.0+.5)); //(int) floor (temp*100.0+.5);

  //now read the poly c hairpin penalty intercept:
  infile>>token;
  while(token != "-->")
    infile>>token;
  infile >> temp;
  cint_ = static_cast<int> ( floor (atof(temp.c_str())*100.0+.5)); //(int) floor (temp*100.0+.5);

  //now read the poly c penalty for a loop of 3:
  infile >>token;
  while(token !=  "-->")
    infile>>token;
  infile >> temp;
  c3_ = static_cast<int> ( floor (atof(temp.c_str())*100.0+.5)); //(int) floor (temp*100.0+.5);

  // Intermolecular initiation free energy 

  infile >>token;
  while(token != "-->")
    infile>>token;
  infile >> temp;
  init_ =  static_cast<int> ( floor (atof(temp.c_str())*100.0+.5)); //(int) floor (temp*100.0+.5);
  
  //now read the GAIL rule indicator
  infile >> token;
  while(token != "-->")
    infile>>token;
  infile>> temp;
  gail_ =  static_cast<int> ( floor (atof(temp.c_str()) +.5)); //(int) floor (temp+.5);
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
 * Input file "tloop.dat". The tetraloop bonus table. 
 * For hairpin loops with four unpaired nucleotides, 
 * this table is consulted. 
 * If the sequence (starting with the 5' last paired
 * nucleotide and finishing with its 3' paired nucleotide) 
 * appears in the table, 
 * the bonus is applied to the hairpin's stability.
 * The file has 30 entries of the type
 * <pre>
 * Seq    Energy 
 ------------- 
 GGGGAC -3.0  
 GGUGAC -3.0  
 ...
 </pre>
 * Each of the six-mer sequences is transformed into an integer
 * by the function tonumi, and entries such as -3.0 are stored
 * as -300 (integers).
*/
void Datatable::input_tloop_dat(const std::string &path) {
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
    infile >> temp;
    tloop_[numoftloops_][1] = (int) floor (100.0*atof(temp)+0.5);
    //std::cout<< "ntloops[0]="<< numoftloops_ << ": "<<base <<  ":"<< tloop_[numoftloops_][0] << std::endl;
    //std::cout<< "ntloops[1]="<< numoftloops_ << ": "<<  tloop_[numoftloops_][1] << std::endl;
    infile >> token;
  }

}

/** read info from dangle.dat 
 *
 * \verbatim
          X                      
  ------------------     
   A    C    G    U       
  ------------------      
      5' --> 3'         
         AX             
         A                  
      3' <-- 5'       
  .     .     .     .   
\endverbatim
* This file shows the energies for dangling bases.
* Presents the free energy of dangling ends (5' or 3'). These energies are used in multibranch and exterior loops.
* In the above example, the X is dangling from the 3' end of the upper strand. This is the case for the first
* four rows. In the last four rows, the X is dangling from the 5' end of the lower strand.
* Indexing. We have l={1,2} and i,j,k={1,2,3,4}. The first four rows have l=1, and the second four rows
* have l=2. The index "i" refers to the row within the bloc of four (e.g., l=1,i?=2 is the second overall row
* and l=2,i=3 refers to the seventh overall row). The index j refers to the column, and the index k refers to the
* position within the individual block (e.g., k=1: A, k=2: C, k=3: G, k=4: U).
*/
void Datatable::input_dangle_dat(const std::string &path) {
  std::ifstream infile(path.c_str());
  char base[110];
  char temp[300];
  if (! infile.is_open() ) {
    std::cerr << "[ERROR] could not initialize " << path << " for I/O" << std::endl;
  }
  std::string token;
  int count, i,j,k,l;

  //add to dangle the case where X (represented as 0) is looked up
  for (l = 1; l <=2; l++){
    for (i = 0; i <=5; i++){
      if ((i!=0)&&(i!=5))
	for (count=1; count <=60; count++)
	  infile >> token;
      for (j=0; j<=5; j++) {
	for (k=0; k<=5; k++) {
	  if ((i==0)||(j==0)||(k==0)) {
	    dangle_[i][j][k][l] = 0;
	  } else if ((i==5)||(j==5)) {
	    dangle_[i][j][k][l] = s_infinity;
	  } else if ((i!=5)&&(j!=5)&&(k==5)) { /* note there was only one "&" (j!=5)&(k==5) in original code! */
	    dangle_[i][j][k][l] = 0;
	  } else {
	    infile >> token;
	    //cout << lineoftext<<"\n";
	    if (token !=  "."){
	      dangle_[i][j][k][l] = static_cast<int> ( floor (100.0*(atof(token.c_str()))+.5));
	    } else {
	      dangle_[i][j][k][l] = s_infinity;
	    }
	  }
	}
      }
    }
  }
}


/**
 * int22.dat
 * //Read the 2x2 internal loops
  //key iloop22[a][b][c][d][j][l][k][m] =
  \verbatim   a j l b
c k m d\endverbatim
  * For instance, for the following:
  * \verbatim
5' ------> 3' 
 U \/ \_/ A 
 G /\  |  C 
3' <------ 5' 
\endverbatim
A=U, B=A C=G D=C. These are the base-paired nucleotides. The inner loop sequence j,l and k,m is analogous.
  */
void Datatable::input_int22_dat(const std::string &path) {
  std::ifstream infile(path.c_str());
  char base[110];
  char temp[300];
  if (! infile.is_open() ) {
    std::cerr << "[ERROR] could not initialize " << path << " for I/O" << std::endl;
  }
  std::string token;
  int count, i,j,k,l,m,a,b,c,d;
  
  for (count=1; count<=340; count++)
    infile >> token; //get past text in file

  for (i=1; i<=36; i++) {//read each of 21 tables
    for (j=1; j<=39; j++) {
      infile >> token; //get past text in file
    }
    strcpy(base, token.c_str());
    strcpy(base+1, "\0");
    a = tonumi(base);
    for (j=1; j<=3; j++)
      infile >> token;
    strcpy(base, token.c_str());
    strcpy(base+1, "\0");
    b = tonumi(base);
    infile>>token;
    strcpy(base, token.c_str());
    strcpy(base+1, "\0");
    c = tonumi(base);
    for (j=1; j<=3; j++)
      infile>> token;
    strcpy(base, token.c_str());
    strcpy(base+1, "\0");
    d = tonumi(base);
    for (j=1; j<=3; j++)
      infile>> token; //get past text in file
    for (j=1; j<=4; j++) {
      for (k=1; k<=4; k++) {
	for (l=1; l<=4; l++) {
	  for (m=1; m<=4; m++) {
	    infile >> temp;
	    iloop22_[a][b][c][d][j][l][k][m] = static_cast<int> (floor(100.0*atof(temp)+0.5));
	  }
	}
      }
    }
  }
}

/**
 * Read from int21.dat,  the 2x1 internal loop data
 */
void Datatable::input_int21_dat(const std::string &path){
  std::ifstream infile(path.c_str());
  char base[110];
  char temp[300];
  if (! infile.is_open() ) {
    std::cerr << "[ERROR] could not initialize " << path << " for I/O" << std::endl;
  }
  std::string token;
  int count, i,j,a,b,c,d,e,f,g;

  for (i=1; i<=58; i++)
    infile >> token; //get past text at top of file
  for (i=1; i<=6; i++) { //read each row of tables
    for (e=1; e<=4; e++) {
      for (j=1; j<=66; j++)
	infile >> token; //get past text in file
      infile >> token;
      strcpy(base, token.c_str());
      strcpy(base+1, "\0");
      a = tonumi(base);
      for (j=1; j<=11; j++)
	infile >> token; //get past text in file
      infile >> token;
      strcpy(base, token.c_str());
      strcpy(base+1, "\0");
      b = tonumi(base);
      for (j=1; j<=35; j++)
	infile >> token; //get past text in file
      for (c=1; c<=4; c++) {
	for (j=1; j<=6; j++) {
	  switch (j) {
	  case 1:
	    f = 1;
	    g = 4;
	    break;
	  case 2:
	    f = 2;
	    g = 3;
	    break;
	  case 3:
	    f = 3;
	    g = 2;
	    break;
	  case 4:
	    f = 4;
	    g = 1;
	    break;
	  case 5:
	    f = 3;
	    g = 4;
	    break;
	  case 6:
	    f = 4;
	    g = 3;
	    break;
	  }
	  for (d=1; d<=4; d++) {
	    infile >> temp;
	    iloop21_[a][b][c][d][e][f][g]= static_cast<int> ( floor(100.0*atof(temp)+0.5));
	  }
	}
      }
    }

  }
}


/**
 *  Read info from coaxial.dat. 
 * add to the stack table the case where X (represented as 0) is looked up:
 * data arrangement of coax: data->coax[a][b][c][d]
  \verbatim 5'bc3'
3'ad
\endverbatim
* The data is arranged like this
\verbatim        
 Y                       Y                       Y                       Y 
  ------------------      ------------------      ------------------      ------------------   
   A    C    G    U        A    C    G    U        A    C    G    U        A    C    G    U  
  ------------------      ------------------      ------------------      ------------------   
      5' --> 3'               5' --> 3'               5' --> 3'               5' --> 3'     
         UX                      UX                      UX                      UX 
         AY                      CY                      GY                      UY 
      3' <-- 5'               3' <-- 5'               3' <-- 5'               3' <-- 5'   
   .     .     .   -1.3    .     .     .     .     .     .     .   -1.0    .     .     .     .   
   .     .   -2.4    .     .     .     .     .     .     .   -1.5    .     .     .     .     .    
   .   -2.1    .   -1.0    .     .     .     .     .   -1.4    .    0.3    .     .     .     .   
 -0.9    .   -1.3    .     .     .     .     .   -0.6    .   -0.5    .     .     .     .     .    
\endverbatim
* It is placed in coax_[j][i][k][l].
* There are four blocks (indexed by j). The overall column is indexed by i. The individual 4x4 block is
* indexed by k and l.
*/
void Datatable::input_coaxial_dat(const std::string &path){
  std::ifstream infile(path.c_str());
  char base[110];
  char temp[300];
  if (! infile.is_open() ) {
    std::cerr << "[ERROR] could not initialize " << path << " for I/O" << std::endl;
  }
  std::string token;
  int count, i,j,k,l; 

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
	    coax_[j][i][k][l]=0;
	  } else if ((i==5)||(j==5)||(k==5)||(l==5)) {
	    coax_[j][i][k][l] = s_infinity;
	  } else {
	    infile >> token;
	    if (token != "."){
	      coax_[j][i][k][l] = static_cast<int> ( floor(100.0*(atof(token.c_str()))+.5));
	      //std::cout << "i=" << i << " j=" << j << " k=" << k << " l="<< l << " coax=" << coax_[j][i][k][l] << std::endl;
	    } else {
	      coax_[j][i][k][l] = s_infinity;
	    }
	  }
	}
      }
    }
  }
}

/**
 * triloop.dat
 * A lookup table for hairpin loops of three. It is applied analogously to the tetraloop table in Tloop.dat.
 * File structure
 *\verbatim
 Seq    Energy
 -------------
CAACG	6.8
GUUAC	6.9
\endverbatim
*/
void Datatable::input_triloop_dat(const std::string &path) {
  std::ifstream infile(path.c_str());
  char base[110];
  char temp[300];
  if (! infile.is_open() ) {
    std::cerr << "[ERROR] could not initialize " << path << " for I/O" << std::endl;
  }
  std::string token;
  int count, i,j,k,l; 

  for (count=1; count<=3; count++)
    infile >> token; //get past text in file
  numoftriloops_=0;
  infile >> token;

  for (count=1; count<=s_maxtloop && !infile.eof(); count++){
    std::cout << token<<"\n";

    numoftriloops_++;
    //strcpy(base, token.c_str());
    base[0]=token[0];
    strcpy(base+1, "\0");
    triloop_[numoftriloops_][0] = tonumi(base);
    //strcpy(base, lineoftext+1);
    base[0]=token[1];
    strcpy(base+1, "\0");
    triloop_[numoftriloops_][0] = triloop_[numoftriloops_][0]+  5*tonumi(base);
    //strcpy(base, lineoftext+2);
    base[0]=token[2];
    strcpy(base+1, "\0");
    triloop_[numoftriloops_][0] = triloop_[numoftriloops_][0]+ 25*tonumi(base);
    //strcpy(base, lineoftext+3);
    base[0]=token[3];
    strcpy(base+1, "\0");
    triloop_[numoftriloops_][0] = triloop_[numoftriloops_][0]+ 125*tonumi(base);
    //strcpy(base, lineoftext+4);
    base[0]=token[4];
    strcpy(base+1, "\0");
    triloop_[numoftriloops_][0] = triloop_[numoftriloops_][0]+ 625*tonumi(base);
    infile >> temp;
    triloop_[numoftriloops_][1] = static_cast<int> (floor (100.0*atof(temp)+0.5));

    //cout << data->triloop[data->numoftriloops][1]<< "  "<<data->triloop[data->numoftriloops][0]<<"\n";

    infile >> token;
  }
}

/**
 * Read info from tstackcoax
 * add to the tstackh table the case where X (represented as 0) is looked up:
 */
void Datatable::input_tstackcoax_dat(const std::string &path) {
  std::ifstream infile(path.c_str());
  char base[110];
  char temp[300];
  if (! infile.is_open() ) {
    std::cerr << "[ERROR] could not initialize " << path << " for I/O" << std::endl;
  }
  std::string token;
  int count, i,j,k,l; 
  
  for (count=1; count<=46; count++)
    infile >> token; //get past text in file
  for (i=0; i<=5; i++) {
    if (!(i==0||i==5))
      for (count=1; count<=60; count++)
	infile >> token;
    for (k=0; k<=5; k++) {
      for (j=0; j<=5; j++) {
	for (l=0; l<=5; l++) {
	  if ((i==0)||(j==0)||(k==0)||(l==0)||(i==5)||(j==5)||(k==5)||(l==5)) {
	    tstackcoax_[i][j][k][l]=0;
	  } else {
	    infile >> token;
	    if (token !=  "."){
	      tstackcoax_[i][j][k][l] = static_cast<int> (floor(100.0*(atof(token.c_str()))+.5));
	      //std::cout << "i="<<i << " j=" << j << " k="<<k<< " l="<<l << " tstackcoax=" << tstackcoax_[i][j][k][l] << std::endl;
	    }  else {
	      tstackcoax_[i][j][k][l] = s_infinity;
	    }
	  }
	}
      }
    }
  }
}

/**
 * The first of the tables for coaxial stacking with an intervening mismatch. This is the stack with the open backbone.
 * add to the tstackh table the case where X (represented as 0) is looked up:
 *\verbatim
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
\todo Do not understand the indexing scheme.
 */
void Datatable::input_coaxstack_dat(const std::string &path) {
  std::ifstream infile(path.c_str());
  char base[110];
  char temp[300];
  if (! infile.is_open() ) {
    std::cerr << "[ERROR] could not initialize " << path << " for I/O" << std::endl;
  }
  std::string token;
  int count, i,j,k,l; 

  for (count=1; count<=46; count++)
    infile >> token; //get past text in file
  for (i=0; i<=4; i++) {
    if (!(i==0||i==5))
      for (count=1; count<=60; count++)
	infile >> token;
    for (k=0; k<=4; k++) {
      for (j=0; j<=4; j++) {
	for (l=0; l<=4; l++) {
	  if ((i==0)||(j==0)||(k==0)||(l==0)||(i==5)||(j==5)||(k==5)||(l==5)) {
	    coaxstack_[i][j][k][l]=0;
	  } else {
	    infile >> token;
	    if (token !=  "."){
	      coaxstack_[i][j][k][l] = static_cast<int> ( floor(100.0*(atof(token.c_str()))+.5));
	      //std::cout << "i="<<i << " j=" << j << " k="<<k<< " l="<<l << " coaxstack=" << coaxstack_[i][j][k][l] << std::endl;
	    } else {
	      coaxstack_[i][j][k][l] = s_infinity;
	    }
	  }
	}
      }
    }
  }
}

/**
 * this function calculates whether a terminal pair i, j requires the end penalty
 */
int Datatable::penalty2(int i, int j) const {
  if (i==4||j==4)
    return auend_;
  else return 0; //no end penalty
}

/**
 * For terminal mismatch stacking in exterior loops, i.e. loops that contain the ends of the sequence.
 *  Read info from tstack 
 * this is the terminal mismatch data used in intermolecular folding
 * add to the tstack table the case where X (represented as 0) is looked up.
 * also add the case where 5 (the intermolecular linker) is looked up,
 * this is actually a dangling end, not a terminal mismatch.
 */
void Datatable::input_tstack_dat(const std::string &path) {
  std::ifstream infile(path.c_str());
  char base[110];
  char temp[300];
  if (! infile.is_open() ) {
    std::cerr << "[ERROR] could not initialize " << path << " for I/O" << std::endl;
  }
  std::string token;
  int count, i,j,k,l;
  
  for (count=1; count<=46; count++)
    infile>> token; //get past text in file
  for (i=0; i<=5; i++) {
    if ((i!=0)&&(i!=5))
      for (count=1; count<=60; count++)
	infile >> token;
    for (k=0; k<=5; k++) {
      for (j=0; j<=5; j++) {
	for (l=0; l<=5; l++) {
	  if ((i==0)||(j==0)||(k==0)||(l==0)) {
	    tstack_[i][j][k][l]=0;
	  } else if ((i==5)||(j==5)) {
	    tstack_[i][j][k][l] = s_infinity;
	  } else if ((k==5)||(l==5)) {
	    //include "5", linker for intermolecular for case of flush ends
	    if ((k==5)&&(l==5)) {//flush end
	      tstack_[i][j][k][l]=0;
	    } else if (k==5) {//5' dangling end
	      //look up number for dangling end
	      tstack_[i][j][k][l] = dangle_[i][j][l][2]+penalty2(i, j);
	    } else if (l==5) {//3' dangling end
	      tstack_[i][j][k][l] = dangle_[i][j][k][1]+penalty2(i, j);
	    }
	  } else {
	    infile>> token;
	    if (token != "."){
	      tstack_[i][j][k][l] =(int) floor (100.0*(atof(token.c_str()))+.5);
	    } else {
	      tstack_[i][j][k][l] = s_infinity;
	    }
	  }
	}
      }
    }
  }
}

/* Read info from tstackm. Used for terminal stacking in a multibranch loop.
 * add to the tstackm table the case where X (represented as 0) is looked up:
 */
void Datatable::input_tstackm_dat(const std::string &path) {
  std::ifstream infile(path.c_str());
  char base[110];
  char temp[300];
  if (! infile.is_open() ) {
    std::cerr << "[ERROR] could not initialize " << path << " for I/O" << std::endl;
  }
  std::string token;
  int count, i,j,k,l;

  for (count=1; count<=46; count++)
    infile >> token; //get past text in file
  for (i=0; i<=4; i++) {
    if (i!=0)
      for (count=1; count<=60; count++)
	infile >> token;
    for (k=0; k<=4; k++) {
      for (j=0; j<=4; j++) {
	for (l=0; l<=4; l++) {
	  if ((i==0)||(j==0)||(k==0)||(l==0)) {
	    tstkm_[i][j][k][l]=0;
	  } else {
	    infile>> token;
	    if (token !=  "."){
	      tstkm_[i][j][k][l] = static_cast<int> (floor(100.0*(atof(token.c_str()))+.5));
	    } else {
	      tstkm_[i][j][k][l] = s_infinity;
	    }
	  }
	}
      }
    }
  }
}

/**
 * Read the 1x1 internal loop data
 * encode the data like:  abc
 *                       def where b-e is a mismatch
 */
void Datatable::input_int11_dat(const std::string &path){
  std::ifstream infile(path.c_str());
  char base[110];
  char temp[300];
  if (! infile.is_open() ) {
    std::cerr << "[ERROR] could not initialize " << path << " for I/O" << std::endl;
  }
  std::string token;
  int count, a,b,c,d,e,f,i,j;

  for (i=1; i<=58; i++)
    infile >> token; //get past text at top of file
  for (i=1; i<=6; i++) { //read each row of table
    if (i==1) {
      a = 1;
      d = 4;
    }
    else if (i==2) {
      a = 2;
      d = 3;
    }
    else if (i==3) {
      a = 3;
      d = 2;
    }
    else if (i==4) {
      a = 4;
      d = 1;
    }
    else if (i==5) {
      a = 3;
      d = 4;
    }
    else {
      a = 4;
      d = 3;
    }
    for (j=1; j<=114; j++)
      infile >> token; //get past text
    for (b=1; b<=4; b++) {
      for (j=1; j<=6; j++) {
      	if (j==1) {
	  c = 1;
	  f = 4;
	}
	else if (j==2) {
	  c = 2;
	  f = 3;
	}
	else if (j==3) {
	  c = 3;
	  f = 2;
	}
	else if (j==4) {
	  c = 4;
	  f = 1;
	}
	else if (j==5) {
	  c = 3;
	  f = 4;
	}
	else {
	  c = 4;
	  f = 3;
	}
	for (e=1; e<=4; e++) {
	  infile >> temp;
	  iloop11_[a][b][c][d][e][f]= static_cast<int>(floor(100.0*atof(temp)+0.5));
	}
      }
    }
  }
}

/**
 * This function is intended for testing only and can be removed at a later time
 * if desired. It takes a 6mer sequence (of the type found in tloop.dat) and
 * looks up the tetraloop energy. It first needs to find the entry that
 * corresponds to the 6mer sequence (there are only 30 of them in the file).
 * If the sequence is not found, return 0. \Note Still trying to understand
 * tloop.dat.
 */
int Datatable::get_tetraloop_energy(const char *seq) const {
  int y=0;
  char base[2];
  base[1]='\0';
  size_t len = strlen(seq);
  if (len != 6)
    return 0;
  base[0]=seq[0];
  y = tonumi(base);
  base[0]=seq[1];
  y = y + 5*tonumi(base);
  base[0]=seq[2];
  y = y + 25*tonumi(base);
  base[0]=seq[3];
  y = y + 125*tonumi(base);
  base[0]=seq[4];
  y = y + 625*tonumi(base);
  base[0]=seq[5];
  y = y + 3125*tonumi(base);
  for (int i=1;i<s_maxtloop;++i) {
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

float Datatable::get_prelog() const {
  return prelog_;
}

int Datatable::get_maxpen() const {
  return maxpen_;
}

int Datatable::get_poppen(unsigned int i) const{
  return poppen_[i];
}

/**
 * This value corresponds to multibranched loops 
 offset in the miscloop file (3.4, should be 340). */
int Datatable::get_constant_multiloop_penalty() const {
  return eparam_[5];
}

/** constant multi-loop penalty for efn2 */
int Datatable::get_constant_efn2_multiloop_penalty() const {
  return efn2a_;
}

/** @return terminal AU penalty (auend_, from miscloop.dat)  */
int Datatable::get_terminal_AU_penalty() const {
  return auend_;
}

int Datatable::get_GU_bonus() const {
  return gubonus_;
}

int Datatable::get_c_hairpin_intercept() const {
  return  cint_;
}
 
int Datatable::get_c_hairpin_slope() const {
  return cslope_;
}
int Datatable::get_c_hairpin_of_3() const{
  return c3_;
}

int Datatable::get_intermolecular_initiation_free_energy() const {
  return init_;
}
  
int Datatable::get_GAIL() const {
  return gail_;
}

/**
 * Return the dangle free energy (dangle.dat).
 * \param i Refers to the row within one of the two blocks of four rows
 * \param j Refers to the column of the blocks
 * \param k position within the individual block (e.g., k=1: A, k=2: C, k=3: G, k=4: U).
 * \param l Refers to the overall four-block (l=1 first four, l=2 second four).
 */
 int Datatable::get_dangle_energy(int i, int j, int k, int l) const {
  return dangle_[i][j][k][l];
}


int Datatable::get_iloop22(int a, int b, int c, int d, int j, int k, int l, int m) const {
  return iloop22_[a][b][c][d][j][k][l][m];
}

int Datatable::get_coaxial_energy(int i, int j, int k, int l) const {
  return coax_[i][j][k][l];
}

int Datatable::get_tstack_coaxial_energy(int i, int j, int k, int l) const {
    return tstackcoax_[i][j][k][l];
  }
