

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
   enum OptionIndex { FOO };
   const option::Descriptor usage[] = {
     { 'c', "enable-foo",  option::ArgType::NONE, 0 },
     {'d', "dddd", option::ArgType::NONE, 0 },
   };
   //program call: "program -c foo --ddd bar";
   const char* argv[] = {"program", "-c", "foo", "--ddd", "bar"};
   int argc=5;
   option::Parser parser(usage, argc, argv);
   std::string cmd = parser.get_command_string();
   
 
   std::string expected = "program -c foo --ddd bar";
   CHECK_STRINGS_EQUAL(expected,cmd);
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

