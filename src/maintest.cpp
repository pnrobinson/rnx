/*
 * A testing framework for the rnx application. 
 *
 * See the files unittest.[h/cpp]. 
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software, to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to permit
 * persons to whom the Software is furnished to do so, subject to the following
 * conditions:
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * 
 * @author Peter Robinson
 * @version 0.1 (22 Jan 2016)
 */

#include "unittest.h"
#include "Sequence.h"
#include "Nussinov.h"
#include "EnergyFunction2.h"
#include "RNAStructure.h"
#include "optionparser.h"

#include <string>
#include <vector>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

// Files with unit test
#include "unittests/readfastatest.cpp"
#include "unittests/nussinovtest.cpp"
#include "unittests/CTtest.cpp"

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



TEST (optiontest1,Parser) {
  std::vector<option::Descriptor> usage = {
     { 'c', "enable-foo",  option::ArgType::NONE, 0 },
     {'d', "ddd", option::ArgType::NONE, 0 },
   };
   //program call: "program -c foo --ddd bar";
   const char* argv[] = {"program", "-c", "foo", "--ddd", "bar"};
   int argc=5;
   option::Parser parser(usage, argc, argv);
   std::string cmd = parser.get_command_string();
   std::string expected = "program -c foo --ddd bar";
   CHECK_STRINGS_EQUAL(expected,cmd);
   expected = "program";
   std::string name = parser.get_program_name();
   CHECK_STRINGS_EQUAL(expected,name);
   bool hasC = parser.has_option('c');
   CHECK(hasC);
   const std::string val = parser.get_value('c');
   std::string expectedstr = "foo";
   CHECK_STRINGS_EQUAL(expectedstr,val);
   bool hasDDD = parser.has_option("ddd");
   CHECK(hasDDD);
   const std::string val2 = parser.get_value("ddd");
   expectedstr = "bar";
   CHECK_STRINGS_EQUAL(expectedstr,val2);
}


TEST (optiontest2,Parser) {
  std::vector<option::Descriptor> usage = {
     { 'a', "alien",  option::ArgType::NONE, "from another planet" },
     { 'b', "brother", option::ArgType::INTEGER, "fratello" },
     { 'c', "cosmos", option::ArgType::FLOAT, "cosmic rays" },
     { 'd', "dumbo", option::ArgType::STRING, "Description" },
   };
   //program call: "program -c foo --ddd bar";
  const char* argv[] = {"analyser", "-a", "martian", "-b", "hermano", "-c", "pluto", "-d", "elephant", "filename"};
  int argc=10;
  option::Parser parser(usage, argc, argv);
  std::string cmd = parser.get_command_string();
  std::string expect = "analyser -a martian -b hermano -c pluto -d elephant filename";
  CHECK_STRINGS_EQUAL(expect,cmd);
  expect = "analyser";
  std::string name = parser.get_program_name();
  CHECK_STRINGS_EQUAL(expect,name);
  bool hasOpt = parser.has_option('a');
  CHECK(hasOpt);
  hasOpt = parser.has_option('b');
  CHECK(hasOpt);
  hasOpt = parser.has_option('c');
  CHECK(hasOpt);
  hasOpt = parser.has_option('d');
  CHECK(hasOpt);
  hasOpt = parser.has_option('e');
  CHECK((!hasOpt));
  hasOpt = parser.has_option('f');
  CHECK((!hasOpt));
  hasOpt = parser.has_option("alien");
  CHECK(hasOpt);
  hasOpt = parser.has_option("brother");
  CHECK(hasOpt);
  hasOpt = parser.has_option("cosmos");
  CHECK(hasOpt);
  hasOpt = parser.has_option("dumbo");
  CHECK(hasOpt);
  hasOpt = parser.has_option("egregious");
  CHECK((!hasOpt));
  const std::string val = parser.get_value('a');
  std::string expectedstr = "martian";
  CHECK_STRINGS_EQUAL(expectedstr,val);
  const std::string val2 = parser.get_value('b');
  expectedstr = "hermano";
  CHECK_STRINGS_EQUAL(expectedstr,val2);
  const std::string val3 = parser.get_value('c');
  expectedstr = "pluto";
  CHECK_STRINGS_EQUAL(expectedstr,val3);
  const std::string val4 = parser.get_value('d');
  expectedstr = "elephant";
  CHECK_STRINGS_EQUAL(expectedstr,val4);
  int n = parser.nonoption_count();
  CHECK_INTS_EQUAL(1,n);
  const std::string nonopt = parser.get_nonoption(0);
  CHECK_STRINGS_EQUAL("filename", nonopt);
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
/* eof */

