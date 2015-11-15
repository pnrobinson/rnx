

#include "unittest.h"
#include "Sequence.h"

#include <string>
#include <vector>

std::string fasta_path="../testdata/KP942435.fasta";

static inline std::string StringFrom(const std::string& value)
{
  return std::string(value.c_str());
}

TEST( Hello, world ) 
{
  std::string s1("Hello"), s2("Hello"), s3("world");
  CHECK_STRINGS_EQUAL(s1, s2);
  CHECK_STRINGS_EQUAL(s2, s1);
  CHECK(s1 != s3);
}

/**
 * Class used to test FASTA Record.cpp
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
  //std::vector<Record> records;
  //bool res= parseFASTA(fasta_path,records);
  int s = records.size();
  CHECK(s==1);
}


/** check we get the correct number of sequences */
TEST(readfasta1,Record)
{
  std::vector<Record> records;
  bool res= parseFASTA(fasta_path,records);
  int s = records.size();
  CHECK(s==1);
}

/** check we get the correct length of the sequence */
TEST(readfasta2,Record)
{
  std::vector<Record> records;
  bool res= parseFASTA(fasta_path,records);
  Record r = records[0];
  unsigned int sz = r.get_size();
  CHECK(sz==226);
}

/** check we get the correct substring starting
 * at position 3 and with length 4.*/
TEST(readfasta3,Record)
{
  std::vector<Record> records;
  bool res= parseFASTA(fasta_path,records);
  Record r = records[0];
  std::string seq = r.substr(3,4);
  CHECK_STRINGS_EQUAL("AGCT",seq);
}

int main(){   
TestResultStdErr result;
TestRegistry::runAllTests(result);
    return (result.getFailureCount());
}

