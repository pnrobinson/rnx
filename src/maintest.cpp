

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

/** check we get the correct number of sequences */
TEST(readfasta1,Record)
{
  std::vector<Record> records = parseFASTA(fasta_path);
  int s = records.size();
  CHECK(s==1);
}

/** check we get the correct length of the sequence */
TEST(readfasta2,Record)
{
  std::vector<Record> records = parseFASTA(fasta_path);
  Record r = records[0];
  unsigned int sz = r.get_size();
  CHECK(sz==226);
}

int main(){   
TestResultStdErr result;
TestRegistry::runAllTests(result);
    return (result.getFailureCount());
}

