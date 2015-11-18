

#include "unittest.h"
#include "Sequence.h"

#include <string>
#include <vector>

std::string fasta_path="../testdata/KP942435.fasta";



TEST( Hello, world ) 
{
  std::string s1("Hello"), s2("Hello"), s3("world");
  CHECK_STRINGS_EQUAL(s1, s2);
  CHECK_STRINGS_EQUAL(s2, s1);
  CHECK(s1 != s3);
}

/**
 * Class used to test FASTA Record.cpp
 * Note that we rely on the FASTA file called KP942435.fasta
 * being present in the directory ../testdata.
 */
class RecordFixtureSetup : public TestSetup {
public:
  void setup() {
    std::string fasta_path="../testdata/KP942435.fasta";
    bool res= parseFASTA(fasta_path,records);
  }

  void teardown() {
    /* no-op */
  }
    
protected:
  std::vector<Record> records;
};




/** check we get the correct number of sequences */
TESTWITHSETUP(RecordFixture, readfasta1)
{
  int s = records.size();
  CHECK(s==1);
}


/** check we get the correct length of the sequence */
TESTWITHSETUP(RecordFixture,readfasta2)
{
  Record r = records[0];
  unsigned int sz = r.get_size();
  CHECK(sz==226);
}

/** check we get the correct substring starting
 * at position 3 and with length 4.*/
TESTWITHSETUP(RecordFixture,readfasta3)
{
  Record r = records[0];
  std::string seq = r.substr(3,4);
  CHECK_STRINGS_EQUAL("AGCT",seq);
}

/** check we get the correct substring starting
 * at position 3 and with length 4.*/
TESTWITHSETUP(RecordFixture,readfasta4)
{
  Record r = records[0];
  std::string rna = r.get_rna();
  std::string seq=rna.substr(3,4);
  CHECK_STRINGS_EQUAL("AGCU",seq);
}

/** check we get the correct acession number
 * >gi|930685873|gb|KP942435.1| Aeromonas hydrophila strain YFG cytotoxic enterotoxin (act) gene, partial cds.*/
TESTWITHSETUP(RecordFixture,readfasta5)
{
  Record r = records[0];
  std::string acc = r.get_accession_number();
  CHECK_STRINGS_EQUAL("KP942435.1",acc);
}

/** check we get the correct acession number
 * >gi|930685873|gb|KP942435.1| Aeromonas hydrophila strain YFG cytotoxic enterotoxin (act) gene, partial cds.*/
TESTWITHSETUP(RecordFixture,readfasta6)
{
  Record r = records[0];
  unsigned int gi = r.get_gi();
  CHECK(930685873u==gi);
}



int main(){   
TestResultStdErr result;
TestRegistry::runAllTests(result);
    return (result.getFailureCount());
}

