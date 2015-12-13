/**
 * unittest.cpp
 */

#include "unittest.h"

#include <exception>
#include <iostream>
#include <sstream>



Test::Test (const std::string& testName) : name (testName) 
{
  TestRegistry::addTest (this);
}


void Test::run (TestResult& result) 
{
  try {
    setup();
    runTest (result);
  } catch (std::exception e) {
    result.addFailure (Failure (e.what(), name, "", 0));	
  } catch (...) {
    result.addFailure (Failure ("Unhandled exception", name, "", 0));	
  }
  teardown();
  result.testWasRun();
}




void TestRegistry::addTest (Test *test) 
{
	instance ().add (test);
}


void TestRegistry::runAllTests (TestResult& result) 
{
	instance ().run (result);
}


TestRegistry& TestRegistry::instance () {
  static TestRegistry registry;
  return registry;
}


void TestRegistry::add (Test *test) {
  tests.push_back (test);
}


void TestRegistry::run (TestResult& result) {
  result.startTests ();
  std::cout << std::endl;
  for (std::vector<Test *>::iterator it = tests.begin (); it != tests.end (); ++it) {
    (*it)->run (result);
    std::cout << ".";
  }
  std::cout << std::endl << std::endl;
  result.endTests ();
}





TestResult::TestResult() : failureCount (0), testCount(0), secondsElapsed(0)
{
  ::time(&startTime_);
}


void TestResult::testWasRun()
{
  testCount++;
}

void TestResult::startTests () 
{
}


void TestResult::addFailure (const Failure & /*failure*/) 
{
  failureCount++;
}

void TestResult::endTests () 
{
    time_t endTime;
    ::time(&endTime);
    secondsElapsed = endTime - startTime_;
}

void TestResultStdErr::addFailure (const Failure & failure) 
{
    TestResult::addFailure(failure);
    std::cerr << failure << std::endl;
}

void TestResultStdErr::endTests () 
{
    TestResult::endTests();
    std::cerr << testCount << " tests run" << std::endl;
    if (failureCount == 1)
      std::cerr << "****** There was one failure";
    else if (failureCount > 1)
        std::cerr << "****** There were " << failureCount << " failures";
    else
      std::cerr << "There were no test failures";
    std::cerr << " (time: " << secondsElapsed << " s)." << std::endl;
}




