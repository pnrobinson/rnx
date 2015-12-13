#ifndef UNITTEST_H
#define UNITTEST_H

/**
 * This is a unit testing framework based on cppunitlite (http://c2.com/cgi/wiki?CppUnitLite) with minor alterations for the rnx project.
 * All of the relevant files have been combined into one header and one source file.
 */

#include <iostream>
#include <math.h>
#include <string>
#include <time.h>
#include <vector>
#include <cstring>

/**
 * Base class from which the Test class will derive.
 */
class TestSetup {
 public:
  virtual void setup() = 0;
  virtual void teardown() = 0;
};




/**
 * This class encapsulates the information that accrues upon a test error.
 */
class Failure
{
 public:
  Failure (std::string theCondition, std::string theTestName, std::string theFileName, long theLineNumber) 
    : condition (theCondition), testName (theTestName), fileName (theFileName), lineNumber (theLineNumber)
  {
    /* no-op */
  }

  std::string condition;
  std::string testName;
  std::string fileName;
  long lineNumber;
};


inline std::ostream& operator<< (std::ostream& stream, const Failure& failure)
{
  stream 
    << failure.fileName.c_str ()
    << " (l. " << failure.lineNumber << "): "
    << failure.testName 
    << " failed: \"" << failure.condition.c_str () << "\" " 
    << std::endl;
  
  return stream;
}



class TestResult;


/**
 * The main class used for creating a unit test.
 */
class Test : public TestSetup
{
 public:
  Test (const std::string& testName);
  virtual void	run (TestResult& result);
  virtual void	runTest (TestResult& result) = 0;
  
 protected:
  std::string name;
};

// The following MACROs are use to define the units tests in client code


#define TEST(name,classUnderTest)\
  class classUnderTest##name##Test : public Test	\
  {							\
  public:								\
    classUnderTest##name##Test () : Test (#name "_Test") {}		\
      void setup() {};							\
      void teardown() {};						\
      void runTest (TestResult& result_);				\
  } classUnderTest##name##Instance;					\
  void classUnderTest##name##Test::runTest (TestResult& result_)	\
  

#define TESTWITHSETUP(name,classUnderTest)				\
  class classUnderTest##name##Test : public Test, name##Setup		\
  {									\
  public:								\
    classUnderTest##name##Test () : Test (#name "Test") {}		\
      void setup() {name##Setup::setup();}				\
      void teardown() {name##Setup::teardown();}			\
      void runTest (TestResult& result_);				\
  } classUnderTest##name##Instance;					\
  void classUnderTest##name##Test::runTest (TestResult& result_)	\
  


#define CHECK(condition) \
    try { \
    if (!(condition)) \
       result_.addFailure (Failure (#condition, name, __FILE__, __LINE__)); \
    } catch(...) { \
        result_.addFailure (Failure ("Unhandled exception", name, __FILE__, __LINE__)); \
    }

#define CHECK_LONGS_EQUAL(expected,actual)\
{\
    try { \
        long _expected = (expected);\
        long _actual = (actual);\
        if (_expected != _actual) {\
            char message [80];\
            sprintf (message, "expected %ld but was: %ld", _expected, _actual);\
            result_.addFailure (Failure (message, name, __FILE__, __LINE__));\
        }\
    } catch(...) { \
    result_.addFailure (Failure ("Unhandled exception", name, __FILE__, __LINE__)); \
    }\
}


#define CHECK_DOUBLES_EQUAL(expected,actual)\
{\
    try { \
        double _expected = (expected);\
        double _actual = (actual);\
        if (fabs ((_expected)-(_actual)) > 0.001) {\
            char message [80];\
            sprintf (message, "expected %lf but was: %lf", (_expected), (_actual));\
            result_.addFailure (Failure (message, name, __FILE__, __LINE__));\
        }\
    } catch(...) { \
    result_.addFailure (Failure ("Unhandled exception", name, __FILE__, __LINE__)); \
    }\
}


#define CHECK_POINTS_EQUAL(expected,actual)\
{\
    try { \
        Point3d _expected = (expected); \
        Point3d _actual = (actual); \
        if (!_actual.Equals(_expected)) { \
            char message [256];\
            sprintf (message, "expected point (%f, %f, %f) but was: (%f, %f, %f)", \
                 _expected.x, _expected.y, _expected.z,\
                 _actual.x, _actual.y, _actual.z); \
            result_.addFailure (Failure (message, name, __FILE__, __LINE__));\
        }\
    } catch(...) { \
    result_.addFailure (Failure ("Unhandled exception", name, __FILE__, __LINE__)); \
    }\
}


#define CHECK_STRINGS_EQUAL(expected,actual)\
{\
    try { \
        std::string _expected(expected);\
        std::string _actual(actual);\
        if (_expected != _actual) {\
            std::string msg = std::string("expected '") + _expected + \
                              std::string("' but was: '") + _actual + "'"; \
            result_.addFailure (Failure (msg.c_str(), name, __FILE__, __LINE__));\
        }\
    } catch(...) { \
        result_.addFailure (Failure ("Unhandled exception", name, __FILE__, __LINE__)); \
    }\
}


#define CHECK_CSTRINGS_EQUAL(expected,actual)\
  {\
    try { \
      int c = strcmp(expected,actual); \
      if (c!=0) { \
	std::string msg=std::string("expected '") + expected +	\
	  std::string("' but was '") + actual + "'";		\
	result_.addFailure(Failure(msg.c_str(), name, __FILE__,__LINE__));\
      } \
    } catch(...) { \
      result_.addFailure (Failure ("Unhandled exception", name, __FILE__, __LINE__)); \
    }\
}


/**
 * This class encapsultes the results of unit testing 
 */
class TestResult
{
 public:
  TestResult ();
  virtual ~TestResult() {};
  virtual void testWasRun ();
  virtual void startTests ();
  virtual void addFailure (const Failure & failure);
  virtual void endTests ();
  int getFailureCount() const { return failureCount; }
  
 protected:
  int failureCount;
  int testCount;
  time_t startTime_;
  int secondsElapsed;
};


class TestResultStdErr : public TestResult
{
 public:
  virtual void addFailure (const Failure & failure);
  virtual void endTests ();
};


class TestResultDebugOut : public TestResult
{
 public:
  virtual void startTests ();
  virtual void addFailure (const Failure & failure);
  virtual void endTests ();
};


/**
 * TestRegistry is a primitive singleton which collects all of the
 * tests in a system and allows them to be executed.  
 */
class TestRegistry
{
 public:
  static void addTest (Test *test);
  static void runAllTests (TestResult& result);
  
 private:
  static TestRegistry&	instance ();
  void add (Test *test);
  void run (TestResult& result);
  
  std::vector<Test *> tests;
};



#endif
