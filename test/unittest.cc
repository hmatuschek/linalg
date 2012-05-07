#include "unittest.hh"
#include <iostream>
#include <sstream>
#include <limits>
#include <cmath>
#include "cputime.hh"

using namespace UnitTest;



/* ********************************************************************************************* *
 * Implementation of TestFailure
 * ********************************************************************************************* */
TestFailure::TestFailure(const std::string &message) throw()
  : message(message)
{
  // pass...
}

TestFailure::~TestFailure() throw()
{
  // Pass...
}

const char *
TestFailure::what() const throw ()
{
  return this->message.c_str();
}




/* ********************************************************************************************* *
 * Implementation of TestCase
 * ********************************************************************************************* */
void
TestCase::setUp()
{
  // Pass...
}


void
TestCase::tearDown()
{
  // Pass...
}


void
TestCase::assertTrue(bool test, const std::string &file, size_t line)
{
  if (! test)
  {
    std::stringstream str;
    str << "Assert failed in " << file << " in line " << line;
    throw TestFailure(str.str());
  }
}


void
TestCase::assertEqual(bool t, bool e, const std::string &file, size_t line)
{
  if (e != t)
  {
    std::stringstream str;
    str << "Expected: " << e << " but got: " << t
        << " in file "<< file << " in line " << line;
    throw TestFailure(str.str());
  }
}


void
TestCase::assertEqual(char t, char e, const std::string &file, size_t line)
{
  if (e != t)
  {
    std::stringstream str;
    str << "Expected: " << e << " but got: " << t
        << " in file "<< file << " in line " << line;
    throw TestFailure(str.str());
  }
}


void
TestCase::assertEqual(short t, short e, const std::string &file, size_t line)
{
  if (e != t)
  {
    std::stringstream str;
    str << "Expected: " << e << " but got: " << t
        << " in file "<< file << " in line " << line;
    throw TestFailure(str.str());
  }
}


void
TestCase::assertEqual(int t, int e, const std::string &file, size_t line)
{
  if (e != t)
  {
    std::stringstream str;
    str << "Expected: " << e << " but got: " << t
        << " in file "<< file << " in line " << line;
    throw TestFailure(str.str());
  }
}


void
TestCase::assertEqual(long t, long e, const std::string &file, size_t line)
{
  if (e != t)
  {
    std::stringstream str;
    str << "Expected: " << e << " but got: " << t
        << " in file "<< file << " in line " << line;
    throw TestFailure(str.str());
  }
}


void
TestCase::assertEqual(unsigned char t, unsigned char e, const std::string &file, size_t line)
{
  if (e != t)
  {
    std::stringstream str;
    str << "Expected: " << e << " but got: " << t
        << " in file "<< file << " in line " << line;
    throw TestFailure(str.str());
  }
}


void
TestCase::assertEqual(unsigned short t, unsigned short e, const std::string &file, size_t line)
{
  if (e != t)
  {
    std::stringstream str;
    str << "Expected: " << e << " but got: " << t
        << " in file "<< file << " in line " << line;
    throw TestFailure(str.str());
  }
}


void
TestCase::assertEqual(unsigned int t, unsigned int e, const std::string &file, size_t line)
{
  if (e != t)
  {
    std::stringstream str;
    str << "Expected: " << e << " but got: " << t
        << " in file "<< file << " in line " << line;
    throw TestFailure(str.str());
  }
}


void
TestCase::assertEqual(unsigned long t, unsigned long e, const std::string &file, size_t line)
{
  if (e != t)
  {
    std::stringstream str;
    str << "Expected: " << e << " but got: " << t
        << " in file "<< file << " in line " << line;
    throw TestFailure(str.str());
  }
}


void
TestCase::assertEqual(float t, float e, const std::string &file, size_t line)
{
  if (e != t)
  {
    std::stringstream str;
    str << "Expected: " << e << " but got: " << t
        << " in file "<< file << " in line " << line;
    throw TestFailure(str.str());
  }
}


void
TestCase::assertEqual(double t, double e, const std::string &file, size_t line)
{
  if (e != t)
  {
    std::stringstream str;
    str << "Expected: " << e << " but got: " << t
        << " in file "<< file << " in line " << line;
    throw TestFailure(str.str());
  }
}


#ifdef UNITTEST_WITH_GINAC
void
TestCase::assertEqual(const GiNaC::ex &t, const GiNaC::ex &e, const std::string &file, size_t line)
{
  if (!e.is_equal(t))
  {
      std::stringstream str;
      str << "Ginac expressions not equal"
          << " in file "<< file << " in line " << line;
      throw TestFailure(str.str());
  }
}
#endif


void
TestCase::assertNear(float t, float e, const std::string &file, size_t line)
{
  if (std::abs(e-t) > std::numeric_limits<float>::epsilon())
  {
    std::stringstream str;
    str << "Expected: " << e << " but got: " << t
        << " in file "<< file << " in line " << line;
    throw TestFailure(str.str());
  }
}


void
TestCase::assertNear(double t, double e, const std::string &file, size_t line)
{
  if (std::abs(e-t) > 10*std::max(std::abs(e),std::abs(t))*std::numeric_limits<double>::epsilon())
  {
    std::stringstream str;
    str << "Expected: " << e << " but got: " << t
        << " in file "<< file << " in line " << line;
    throw TestFailure(str.str());
  }
}




/* ********************************************************************************************* *
 * Implementation of TestSuite
 * ********************************************************************************************* */
TestSuite::TestSuite(const std::string &desc)
  : description(desc)
{
  // pass...
}


TestSuite::~TestSuite()
{
  // Free callers:
  for (std::list<TestCallerInterface *>::iterator caller=this->tests.begin();
       caller != this->tests.end(); caller++) {
    delete *caller;
  }
}


void
TestSuite::addTest(TestCallerInterface *test)
{
  this->tests.push_back(test);
}


const std::string &
TestSuite::getDescription()
{
  return this->description;
}


TestSuite::iterator
TestSuite::begin()
{
  return this->tests.begin();
}


TestSuite::iterator
TestSuite::end()
{
  return this->tests.end();
}




/* ********************************************************************************************* *
 * Implementation of TestSuite
 * ********************************************************************************************* */
TestRunner::TestRunner(std::ostream &stream)
  : stream(stream)
{
  // Pass...
}


TestRunner::~TestRunner()
{
  // Free suites:
  for (std::list<TestSuite *>::iterator suite = this->suites.begin();
       suite != this->suites.end(); suite++)
  {
    delete *suite;
  }
}


void
TestRunner::addSuite(TestSuite *suite)
{
  this->suites.push_back(suite);
}


void
TestRunner::operator ()()
{
  size_t tests_run = 0;
  size_t tests_failed = 0;
  size_t tests_error = 0;

  for (std::list<TestSuite *>::iterator suite = this->suites.begin();
       suite != this->suites.end(); suite++)
  {
    // Dump Suite description
    this->stream << "Suite: " << (*suite)->getDescription() << std::endl;

    // For each test in suite:
    for (TestSuite::iterator test = (*suite)->begin(); test != (*suite)->end(); test++)
    {
      this->stream << " test: " << (*test)->getDescription() << ": ";

      try
      {
        tests_run++;
        CpuTime clock; clock.start();
        // Run test
        (**test)();
        this->stream << " ok (" << clock.stop() << "s)" << std::endl;
      }
      catch (TestFailure fail)
      {
        this->stream << " fail" << std::endl;
        this->stream << "   reason: " << fail.what() << std::endl;
        tests_failed++;
      }
      catch (std::exception &err)
      {
        this->stream << " exception" << std::endl;
        this->stream << "   what(): " << err.what() << std::endl;
        tests_error++;
      }
    }

    this->stream << std::endl;
  }

  this->stream << "Summary: " << tests_failed << " tests failed out of "
               << tests_run - tests_error
               << " (" << 100. * float((tests_run-tests_failed-tests_error))/(tests_run-tests_error)
               << "% passed)." << std::endl << "         Where "
               << tests_error << " tests produced errors." << std::endl;
}

