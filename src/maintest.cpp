

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
    delete rnastruct;
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

/** check we get the correct label of the RA7680 structure */
TESTWITHSETUP(RNAStructureFixture,intermolecular)
{
  bool intermol = rnastruct->intermolecular();// Note using one-based numbering!
  bool expected = false;
  CHECK(expected==intermol);
}

/** check erg1 (see latex tutorial) */
TESTWITHSETUP(RNAStructureFixture,erg1)
{
  const char *dir = "../dat";
  Datatable dattab(dir);
  int val = dattab.erg1(2,71,3,70,rnastruct);
  int expected = -150;
  CHECK_INTS_EQUAL(expected, val);
}


/** check erg3 (see latex tutorial) */
TESTWITHSETUP(RNAStructureFixture,erg3)
{
  const char *dir = "../dat";
  Datatable dattab(dir);
  int dbl = 1;
  int e3 = dattab.erg3(2,50,rnastruct,dbl);
  // expect "infinity"
  CHECK(e3>9999);
  // Now simulate the situation in which the haipin comrpises the
  // entire sequence 1..n, in this case the energy should be infinity
  int n = rnastruct->get_number_of_bases();
  dbl=42; // not 1, and not 5--not sure what appropriate code is
  e3 = dattab.erg3(1,n,rnastruct,dbl);
  // expect "infinity"
  CHECK(e3>9999);
  // Now check for hairpin loops longer than 30 nt.
  // To simulate this with the sequence of our tRNA
  // take i=10 (G), j=56 (C), then i+1=11 (C), j-1=55(U)
  // the corresponding indices to tstackh are  [G][C][C][U] or [3][2][2][4]
  // Looking this value up in tstackh.dat, I see -0.5 or -50
  // We need to add this to the value of hairpin_[30]
  // This is parsed from loop dat, and can be seen in the "HAIRPIN" column of
  // that file, row 30 to be 7.7 or 770
  // Therefore, we expect -50+770+ 107.857764*log(n/30)
  // where n is the size of the loop, j-i+1=56-10+1=47
  //ln(47/30)=0.4490
  // 107.857764*log(n/30) = 43
  // Therefore, we expect -50+770+43=763
  // energy = tstkh_[ct->numseq(i)][ct->numseq(j)][ct->numseq(i+1)][ct->numseq(j-1)]
  //    + hairpin_[30]+loginc+eparam_[4];
  int expected = 763;
  e3 = dattab.erg3(10,56,rnastruct,dbl);
  CHECK_INTS_EQUAL(expected,e3);
  // Now check for hairpin loops shorter than 3 nucleotides. This is prohibited by
  // nearest neighbor theory and should return infinity.
  // The following means that that bases 10 and 12 are paired and there is a "loop"
  // of size 1 between them.
  e3 = dattab.erg3(10,12,rnastruct,dbl);
  CHECK(e3>9999);
  // The following means that that bases 10 and 13 are paired and there is a "loop"
  // of size 2 between them.
  e3 = dattab.erg3(10,13,rnastruct,dbl);
  CHECK(e3>9999);
  // Simulate a tetraloop over positions 13..18: CAGAUG
  // This should be a simple lookup in tstkh
  // i=13 (C), i+1=14 (A), j=18 (G) j-1=17 (U)
  // indices are [2][3][1][4]
  // The file tstackh.dat shows -1.8 at this
  //position.
  expected = 380;
  e3 = dattab.erg3(13,18,rnastruct,dbl);
  CHECK_INTS_EQUAL(expected,e3);
  // Simulate a triloop over positions   31=CUUGG=35
  // see latex document for details.
  expected = 570;
  e3 = dattab.erg3(31,35,rnastruct,dbl);
  CHECK_INTS_EQUAL(expected,e3);
  // Check a longer loop, from position 13 (C)-14(A)-...21(A)-22(G)
  // index[2][3][1][1]. -1.5 or -150
  // then we have a hairpin penalty value for a length of 22-1311=8 of 5.6 or 560
  // Note that size does not include the two paired bases i and j
  // expected is thus -150+560=410
  expected = 410;
  e3 = dattab.erg3(13,22,rnastruct,dbl);
  CHECK_INTS_EQUAL(expected,e3);
}


/** check erg2 (see latex tutorial) */
TESTWITHSETUP(RNAStructureFixture,erg2)
{
  const char *dir = "../dat";
  Datatable dattab(dir);
  int dbl = 42;
  int e2;
  int a=42;
  int b=42; // a and b have no meaning for the following test.
  int numbases = rnastruct->get_number_of_bases();
  //A loop cannot contain the ends of the sequence (where one of the indices is more than the seqlen - duh
  //Test that function returns "infinite" energy
  e2 =  dattab.erg2(10, numbases+1, 25, 30, rnastruct, a, b);
  CHECK(e2>9999);
  e2 =  dattab.erg2(10, 15, numbases+1, 30, rnastruct, a, b);
  CHECK(e2>9999);
}



TEST (getnumbases,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  std::string ct_file = "../testdata/SNORA17.ct";
  RNAStructure rnastruct(ct_file);
  int numbases = rnastruct.get_number_of_bases();
  int expected = 132;
  CHECK_INTS_EQUAL(expected,numbases);
}


TEST (rnabulge_1,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  std::string ct_file = "../testdata/SNORA17.ct";
  RNAStructure rnastruct(ct_file);
  int dbl = 42;
  int e2;
  int a=42;
  int b=42; // a and b have no meaning for the following test.
  // Note that base $U^{108}$ is a bulge base, and thus
  //there is a bulge loops formed by base pairs $i.j$=($G^{89}$,$C^{109}$) and $i',j'$=($U^{90}$,$A^{107}$)
  // We are thus testing a bulge loop of size 1.
  // Bulge loops up to size 30 are assigned free energies from the loop file
  // Looking up the entry for "BULGE" in loop.dat for size=1, we see 3.8 (i.e., we expect 380)
  // we also add the corresponding entry from stack.dat stack_[G][C][U][A] = stack_[3][2][4][1]=-220
  // thus we expect 380-220=160
  int i=89;
  int j=109;
  int ip=90;
  int jp=107;
  e2 = dattab.erg2(i,j,ip,jp, &rnastruct, a, b);
  int expected = 160;
  CHECK_INTS_EQUAL(expected,e2);
  // Now check a bulge of length 3
  i=15;
  j=47;
  ip = 19;
  jp = 46;
  // The energy for a size-3 bulge is calculate directly from bulge_[3]
  // There is an additional penalty of 50 for the G:U base pairing
  // (auend_).
  e2 = dattab.erg2(i,j,ip,jp, &rnastruct, a, b);
  expected = 370;
  CHECK_INTS_EQUAL(expected,e2);
}

/**
 * Testing inner loop 2x2 with mir283 structure
 */
TEST (iloop22_mir283,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  std::string ct_file = "../testdata/mir283.ct";
  RNAStructure rnastruct(ct_file);
  int dbl = 42;
  int e2;
  int a=42;
  int b=42; // a and b have no meaning for the following test.
  int i=10;
  int j=94;
  int ip=13;
  int jp=91;
  // The energy for the 2x2 innerloop in microrna-282
  // see latex document for the details.
  e2 = dattab.erg2(i,j,ip,jp, &rnastruct, a, b);
  int expected = 120;
  CHECK_INTS_EQUAL(expected,e2);
}


/**
 * Testing the u3 file
 * Note that
 * (1) There are five structures (try $ grep dG u3.ct)
 */
TEST (u3,EnergyFunction2) {
  std::string ct_file = "../testdata/u3.ct";
  RNAStructure rnastruct(ct_file);
  int n = rnastruct.get_number_of_structures();
  int expected = 5;
  CHECK_INTS_EQUAL(expected,n);
  std::string label = rnastruct.get_ith_label(1);
  std::string expected_str = "dG = -82.50 [Initially -85.70] AAPY01491858 1/243-27";
  CHECK_STRINGS_EQUAL(expected_str,label);


}
/**
 * Testing inner loop 2x1 with u3 structure
 */
TEST (iloop21_u3,EnergyFunction2) {
  const char *dir = "../dat";
  Datatable dattab(dir);
  std::string ct_file = "../testdata/u3.ct";
  std::cout << "About to input u3\n";
  RNAStructure rnastruct(ct_file);
  int dbl = 42;
  int e2;
  int a=42;
  int b=42; // a and b have no meaning for the following test.
  int i=19;
  int j=50;
  int ip=22;
  int jp=48;
  // The energy for the 2x2 innerloop in microrna-282
  // see latex document for the details.
  e2 = dattab.erg2(i,j,ip,jp, &rnastruct, a, b);
  int expected = 120;
  CHECK_INTS_EQUAL(expected,e2);
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

