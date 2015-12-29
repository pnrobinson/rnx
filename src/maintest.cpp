

#include "unittest.h"
#include "Sequence.h"
#include "Nussinov.h"
#include "EnergyFunction2.h"
#include "RNAStructure.h"

#include <string>
#include <vector>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

// Files with unit test
#include "unittests/readfastatest.cpp"


/** Handler to print a stack trace if there is a SEGFAULT */
void handler(int sig) {
  void *array[30];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 30);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  fprintf(stderr, "Consider  addr2line 0x435607 -e maintest\n" );
  exit(1);
}



class RNAStructureFixtureSetup : public TestSetup {
  public:
  void setup() {
    std::string ct_file = "../testdata/RA7680.ct";
    rnastruct=new RNAStructure(ct_file);
   
  }

  void teardown() {
    /* nothing to do */
  }
  protected:
    RNAStructure *rnastruct;
};



class NussinovFixtureSetup: public TestSetup {
public:
  void setup() {
    std::string gb_path="../testdata/NM_000518.gb";
    std::vector<Record> records;
    bool res= parseGenBank(gb_path,records);
    Record r = records[0];
    std::string rna = r.get_rna();
    nuss = new Nussinov(rna.c_str());
  }
  void teardown() {
    delete nuss;
  }
protected:
  Nussinov * nuss;
};




/** check we get the correct length of the HBB sequence */
TESTWITHSETUP(NussinovFixture,size1)
{
  unsigned int sz = nuss->get_len();
  CHECK(sz==626);
}



TEST (fold1, Nussinov) {
  const char * rna= "GAAAC";
  Nussinov * nuss = new Nussinov(rna);
  const char * folded = nuss->fold_rna();
  CHECK_CSTRINGS_EQUAL("(...)",folded);
  delete nuss;
}


TEST (fold2, Nussinov) {
  const char * rna= "GAAAAC";
  Nussinov * nuss = new Nussinov(rna);
  const char * folded = nuss->fold_rna();
  CHECK_CSTRINGS_EQUAL("(....)",folded);
  delete nuss;
}

TEST (fold3, Nussinov) {
  const char * rna= "CAAAAAG";
  Nussinov * nuss = new Nussinov(rna);
  const char * folded = nuss->fold_rna();
  CHECK_CSTRINGS_EQUAL("(.....)",folded);
  delete nuss;
}

TEST (fold4, Nussinov) {
  const char * rna= "CGAAUCG";
  Nussinov * nuss = new Nussinov(rna);
  const char * folded = nuss->fold_rna();
  CHECK_CSTRINGS_EQUAL("((...))",folded);
  delete nuss;
}


TEST (fold5, Nussinov) {
  const char * rna= "CGAAACGU";
  Nussinov * nuss = new Nussinov(rna);
  const char * folded = nuss->fold_rna();
  CHECK_CSTRINGS_EQUAL("((...)).",folded);
  delete nuss;
}



TEST (fold6, Nussinov) {
  const char * rna= "GGGAAAUCC";
  Nussinov * nuss = new Nussinov(rna);
  const char * folded = nuss->fold_rna();
  CHECK_CSTRINGS_EQUAL("(((...)))",folded);
  delete nuss;
}



TEST (fold7, Nussinov) {
  const char * rna= "ACUCGAUCCGAG";
  Nussinov * nuss = new Nussinov(rna);
  const char * folded = nuss->fold_rna();
  CHECK_CSTRINGS_EQUAL(".((((...))))",folded);
  delete nuss;
}



TEST (getdatadir,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  CHECK_CSTRINGS_EQUAL(dir,dattab.get_data_dir());
}

/**
 * This tests whether we get the destabilising energies
 * for internal loops. There are 30 values in the table, which
 * comes from loop.dat and looks like this:
 * SIZE         INTERNAL            BULGE            HAIRPIN 
------------------------------------------------------- 
1                .                3.8                 . 
2                .                2.8                 . 
3                .                3.2               5.7   
4              1.7                3.6               5.6  
5              1.8                4.0               5.6  
*/
TEST (getdestabilizingenergy_loop1,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  int d = dattab.get_destabilizing_energy_internal_loop(1);
  CHECK(d>1000); // should be infinity.
  d = dattab.get_destabilizing_energy_internal_loop(2);
  CHECK(d>1000); // should be infinity.
  d = dattab.get_destabilizing_energy_internal_loop(3);
  CHECK(d>1000); // should be infinity.
  d = dattab.get_destabilizing_energy_internal_loop(4);
  CHECK_INTS_EQUAL(170,d);
  d = dattab.get_destabilizing_energy_internal_loop(5);
  CHECK_INTS_EQUAL(180,d);
  d = dattab.get_destabilizing_energy_internal_loop(6);
  CHECK_INTS_EQUAL(200,d);
  d = dattab.get_destabilizing_energy_internal_loop(7);
  CHECK_INTS_EQUAL(220,d);
   d = dattab.get_destabilizing_energy_internal_loop(8);
  CHECK_INTS_EQUAL(230,d);
   d = dattab.get_destabilizing_energy_internal_loop(9);
  CHECK_INTS_EQUAL(240,d);
   d = dattab.get_destabilizing_energy_internal_loop(10);
  CHECK_INTS_EQUAL(250,d);
   d = dattab.get_destabilizing_energy_internal_loop(11);
  CHECK_INTS_EQUAL(260,d);
   d = dattab.get_destabilizing_energy_internal_loop(12);
  CHECK_INTS_EQUAL(270,d);
  // ...
   d = dattab.get_destabilizing_energy_internal_loop(29);
  CHECK_INTS_EQUAL(360,d);
  d = dattab.get_destabilizing_energy_internal_loop(30);
  CHECK_INTS_EQUAL(370,d);
}

/**
 * This tests whether we get the destabilising energies
 * for bulge loops. There are 30 values in the table, which
 * comes from loop.dat and looks like this:
 * SIZE         INTERNAL            BULGE            HAIRPIN 
------------------------------------------------------- 
1                .                3.8                 . 
2                .                2.8                 . 
3                .                3.2               5.7   
4              1.7                3.6               5.6  
5              1.8                4.0               5.6  
...
28             3.6                6.0               7.6  
29             3.6                6.0               7.6  
30             3.7                6.1               7.7  
*/
TEST (getdestabilizingenergy_loop2,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  int d = dattab.get_destabilizing_energy_bulge_loop(1);
  CHECK_INTS_EQUAL(380,d); // 3.8
  d = dattab.get_destabilizing_energy_bulge_loop(2);
  CHECK_INTS_EQUAL(280,d); 
  d = dattab.get_destabilizing_energy_bulge_loop(3);
  CHECK_INTS_EQUAL(320,d);
   d = dattab.get_destabilizing_energy_bulge_loop(28);
  CHECK_INTS_EQUAL(600,d);
   d = dattab.get_destabilizing_energy_bulge_loop(29);
  CHECK_INTS_EQUAL(600,d);
   d = dattab.get_destabilizing_energy_bulge_loop(30);
  CHECK_INTS_EQUAL(610,d); 
}

/**
 * This tests whether we get the destabilising energies
 * for hairpin loops. There are 30 values in the table, which
 * comes from loop.dat and looks like this:
 * SIZE         INTERNAL            BULGE            HAIRPIN 
------------------------------------------------------- 
1                .                3.8                 . 
2                .                2.8                 . 
3                .                3.2               5.7   
4              1.7                3.6               5.6  
5              1.8                4.0               5.6  
...
28             3.6                6.0               7.6  
29             3.6                6.0               7.6  
30             3.7                6.1               7.7  
*/
TEST (getdestabilizingenergy_loop3,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  int d = dattab.get_destabilizing_energy_hairpin_loop(1);
  CHECK(d>9999); // should be infinity 
  d = dattab.get_destabilizing_energy_hairpin_loop(2);
  CHECK(d>9999); // should be infinity 
  d = dattab.get_destabilizing_energy_hairpin_loop(3);
  CHECK_INTS_EQUAL(570,d);
   d = dattab.get_destabilizing_energy_hairpin_loop(28);
  CHECK_INTS_EQUAL(760,d); 
  d = dattab.get_destabilizing_energy_hairpin_loop(29);
  CHECK_INTS_EQUAL(760,d); 
  d = dattab.get_destabilizing_energy_hairpin_loop(30);
  CHECK_INTS_EQUAL(770,d);
}

TEST (getstackenergy_1,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  int d = dattab.get_stack_energy(1,1,1,1);
  CHECK(d>9999); // should be infinity
  d = dattab.get_stack_energy(1,4,1,4);
  CHECK_INTS_EQUAL(-90,d); // AX/UY,row 1 column 4: -0.9
  d = dattab.get_stack_energy(1,4,2,3);
  CHECK_INTS_EQUAL(-220,d); // AX/UY,row 2 column 3: -2.2
  d = dattab.get_stack_energy(4,3,3,2);
  CHECK_INTS_EQUAL(-140,d); // AX/UY,row 2 column 3: -2.2
}

TEST (getstackenergy_2,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  int d = dattab.get_tstackh_energy(1,1,1,1);
  CHECK(d>9999); // should be infinity
  d = dattab.get_tstackh_energy(1,4,1,4);
  CHECK_INTS_EQUAL(-30,d); // AX/UY,row 1 column 4: -0.3
  d = dattab.get_tstackh_energy(1,4,2,3);
  CHECK_INTS_EQUAL(-150,d); // AX/UY,row 2 column 3: -2.2
  d = dattab.get_tstackh_energy(4,3,3,2);
  CHECK_INTS_EQUAL(-120,d); // AX/UY,row 2 column 3: -2.2
}

TEST (getstackenergy_3,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  int d = dattab.get_tstacki_energy(1,1,1,1);
  CHECK(d>9999); // should be infinity
  d = dattab.get_tstacki_energy(1,4,1,4);
  CHECK_INTS_EQUAL(70,d); // AX/UY,row 1 column 4: -0.3
  d = dattab.get_tstacki_energy(1,4,2,3);
  CHECK_INTS_EQUAL(70,d); // AX/UY,row 2 column 3: -2.2
  d = dattab.get_tstacki_energy(4,3,3,2);
  CHECK_INTS_EQUAL(70,d); // AX/UY,row 2 column 3: -2.2
}

/*
 * Need to add more tests, do not understand format yet!
 */
TEST (tloop,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  int d = dattab.get_tetraloop_energy("GGGGAC");
  CHECK_INTS_EQUAL(-300,d); //  -3.0
  d = dattab.get_tetraloop_energy("CGAAAG");
  CHECK_INTS_EQUAL(-300,d); //  -3.0
  d = dattab.get_tetraloop_energy("CGAAGG");
  CHECK_INTS_EQUAL(-250,d); //  -2.5
  d = dattab.get_tetraloop_energy("UGGAAA");
  CHECK_INTS_EQUAL(-150,d); // -1.5
}

/*
 * An item from miscloop.dat
 */
TEST (prelog,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  float pl = dattab.get_prelog();
  float expected = 100.0*1.07857764; /* from miscloop.dat */
  CHECK_DOUBLES_EQUAL(expected,pl);
}

/*
 * An item from miscloop.dat
 */
TEST (maxpen,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  int maxpen = dattab.get_maxpen();
  int expected = 300; /* from miscloop.dat */
  CHECK_INTS_EQUAL(expected,maxpen);
}

/*
 * An item from miscloop.dat
 */
TEST (poppen,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  int poppen = dattab.get_poppen(1);
  int expected = 50; /* from miscloop.dat */
  CHECK_INTS_EQUAL(expected,poppen);
  poppen = dattab.get_poppen(2);
  CHECK_INTS_EQUAL(expected,poppen);
  poppen = dattab.get_poppen(3);
  CHECK_INTS_EQUAL(expected,poppen);
  poppen = dattab.get_poppen(4);
  CHECK_INTS_EQUAL(expected,poppen);
}

TEST(get_constant_multiloop_penalty,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  int mlp = dattab.get_constant_multiloop_penalty();
  int expected = 340; /* from miscloop.dat */
  CHECK_INTS_EQUAL(expected,mlp);
  mlp = dattab.get_constant_efn2_multiloop_penalty();
  expected = 1010; /* from miscloop.dat */
  CHECK_INTS_EQUAL(expected,mlp);
}

TEST(auend,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  int tp = dattab.get_terminal_AU_penalty();
  int expected = 50; /* from miscloop.dat */
  CHECK_INTS_EQUAL(expected,tp);
}

TEST(gubonus,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  int gub = dattab.get_GU_bonus();
  int expected = -220; /* from miscloop.dat */
  CHECK_INTS_EQUAL(expected,gub);
}


TEST(c_hairpin,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  int cint = dattab.get_c_hairpin_intercept();
  int expected = 160; // 1.6
  CHECK_INTS_EQUAL(expected,cint);
  int cslope = dattab.get_c_hairpin_slope();
  expected = 30; // 30
  CHECK_INTS_EQUAL(expected,cslope);
  int c3 = dattab.get_c_hairpin_of_3();
  expected = 140; //1.4
  CHECK_INTS_EQUAL(expected,c3);
}


TEST(get_intermolecular_initiation_free_energy,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  int iife = dattab.get_intermolecular_initiation_free_energy();
  int expected = 410; // 4.1
  CHECK_INTS_EQUAL(expected,iife);
}

/* still miscloop.dat) */
TEST(gail,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  int gail = dattab.get_GAIL();
  int expected = 1; // 1= on, 0=off
  CHECK_INTS_EQUAL(expected,gail);
}


/* dangle.dat). 
* We have l={1,2} and i,j,k={1,2,3,4}. The first four rows have l=1, and the second four rows
* have l=2. The index "i" refers to the row within the bloc of four (e.g., l=1,i?=2 is the second overall row
* and l=2,i=3 refers to the seventh overall row). The index j refers to the column, and the index k refers to the
* position within the individual block (e.g., k=1: A, k=2: C, k=3: G, k=4: U).
* dangle_[i][j][k][l] 
*/
TEST(dangle,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  int d = dattab.get_dangle_energy(1,1,1,1);
  CHECK(d>9999); // should be infinity
  d = dattab.get_dangle_energy(1,4,4,1);
  int expected = -60; // -0.6
  CHECK_INTS_EQUAL(expected,d);
  d = dattab.get_dangle_energy(3,2,2,1);
  expected = -40; // -0.4
  CHECK_INTS_EQUAL(expected,d);
  d = dattab.get_dangle_energy(4,3,1, 2);
  expected = -30; // -0.3
  CHECK_INTS_EQUAL(expected,d);
}

TEST(iloop2,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  int a = 1; // "A"
  int b = 1; // "A"
  int c = 4; // "U"
  int d = 4; // "U"
  int il22 = dattab.get_iloop22(a,b,c,d, 1,1,1,1);
  int expected = 280; // 2.9
  CHECK_INTS_EQUAL(expected,il22);
  il22 = dattab.get_iloop22(a,b,c,d, 1,1,1,2);
  expected = 230; // 2.3
  CHECK_INTS_EQUAL(expected,il22);
  il22 = dattab.get_iloop22(a,b,c,d, 1,1,1,3);
  expected = 170; // 1.7
  CHECK_INTS_EQUAL(expected,il22);

}

/**
 * See input_coax_dat for arrangement of indices.
 */
TEST(coaxial,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  
  int ce = dattab.get_coaxial_energy(3, 2, 4,1);
  int expected = -210; // -2.1, 
  CHECK_INTS_EQUAL(expected,ce);
  ce = dattab.get_coaxial_energy(4, 3, 4,1);
  expected = -140; // -1.4, 
  CHECK_INTS_EQUAL(expected,ce);
}

/**
 * See tstackcoax_dat for arrangement of indices.
 * Do not understand indexing, come back to this!
 */
TEST(tstack_coaxial,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  //i=1 j=1 k=2 l=1 tstackcoax=-110
  int ce = dattab.get_tstack_coaxial_energy(1, 1, 2,1);
  int expected = -110; // -1.1, 
  CHECK_INTS_EQUAL(expected,ce);
}





/** check we get the correct length of the RA7680 sequence */
TESTWITHSETUP(RNAStructureFixture,basics)
{
  int nbase = rnastruct->get_number_of_bases();
  int expected = 76;
  CHECK_INTS_EQUAL(expected, nbase);
  int n = rnastruct->get_number_of_structures();
  expected = 1;
  CHECK_INTS_EQUAL(expected,n);
}



/** check we get the correct label of the RA7680 structure */
TESTWITHSETUP(RNAStructureFixture,labels)
{
  std::string lbl = rnastruct->get_ith_label(1);// Note using one-based numbering!
  std::string expected = "ENERGY = 0.1  RA7680";
  CHECK_STRINGS_EQUAL(expected, lbl);
}



/** check we get the correct dot-paren of the RA7680 CT file.
 * This is what I get from RNAstructure ct2dot Results
 * at http://rna.urmc.rochester.edu/RNAstructureWeb
 * >   ENERGY = 0.1  RA7680
GGGGGCGUAGCUCAGAUGGUAGAGCGCUCGCUUGGCGUGUGAGAGGUACCGGGAUCGAUACCCGGCGCCUCCACCA
(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))....
*/
TESTWITHSETUP(RNAStructureFixture,dotparen)
{
  std::string dpar = rnastruct->get_dot_parens_structure(1);// Note using one-based numbering!
  std::string expected = "(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))....";
  CHECK_STRINGS_EQUAL(expected, dpar);
}



 
int main(){   
  TestResultStdErr result;
  signal(SIGSEGV, handler);
  signal(SIGABRT, handler);
  try {
    TestRegistry::runAllTests(result);
    return (result.getFailureCount());
  } catch(const std::exception &e) {
    std::cout << "Caught exception: " << e.what() << std::endl;
  }

}

