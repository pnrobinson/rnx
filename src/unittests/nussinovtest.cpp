

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
