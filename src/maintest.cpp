

#include "unittest.h"
#include "Sequence.h"
#include "Nussinov.h"
#include "EnergyFunction2.h"

#include <string>
#include <vector>

// Files with unit test
#include "unittests/readfastatest.cpp"







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



int main(){   
  TestResultStdErr result;
  TestRegistry::runAllTests(result);
  return (result.getFailureCount());
}

