/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __FLUC_UNITTEST_HH__
#define __FLUC_UNITTEST_HH__

#include <list>
#include <string>

#ifdef UNITTEST_WITH_GINAC
#include <ginac/ginac.h>
#endif

namespace UnitTest {

class TestFailure : public std::exception
{
protected:
  std::string message;

public:
  TestFailure(const std::string &message) throw();
  virtual ~TestFailure() throw();

  const char *what() const throw();
};



class TestCase
{
public:
  virtual void setUp();
  virtual void tearDown();

  void assertTrue(bool test, const std::string &file, size_t line);
  void assertEqual(bool t, bool e, const std::string &file, size_t line);
  void assertEqual(char t, char e, const std::string &file, size_t line);
  void assertEqual(short t, short e, const std::string &file, size_t line);
  void assertEqual(int t, int e, const std::string &file, size_t line);
  void assertEqual(long t, long e, const std::string &file, size_t line);
  void assertEqual(unsigned char t, unsigned char e, const std::string &file, size_t line);
  void assertEqual(unsigned short t, unsigned short e, const std::string &file, size_t line);
  void assertEqual(unsigned int t, unsigned int e, const std::string &file, size_t line);
  void assertEqual(unsigned long t, unsigned long e, const std::string &file, size_t line);
  void assertEqual(float t, float e, const std::string &file, size_t line);

#ifdef UNITTEST_WITH_GINAC
  void assertEqual(const GiNaC::ex &t,const GiNaC::ex &e, const std::string &file, size_t line);
#endif

  void assertEqual(double t, double e, const std::string &file, size_t line);
  void assertNear(float t, float e, const std::string &file, size_t line);
  void assertNear(double t, double e, const std::string &file, size_t line);
};



class TestCallerInterface
{
protected:
  std::string description;

public:
  TestCallerInterface(const std::string &desc)
    : description(desc)
  {
    // Pass...
  }

  virtual const std::string &getDescription()
  {
    return this->description;
  }

  virtual void operator() () = 0;
};


template <class T>
class TestCaller : public TestCallerInterface
{
protected:
  void (T::*function)(void);

public:
  TestCaller(const std::string &desc, void (T::*func)(void))
    : TestCallerInterface(desc), function(func)
  {
    // Pass...
  }

  virtual void operator() ()
  {
    // Create new test:
    T *instance = new T();

    // Call test
    instance->setUp();
    (instance->*function)();
    instance->tearDown();

    // free instance:
    delete instance;
  }
};


class TestSuite
{
public:
  typedef std::list<TestCallerInterface *>::iterator iterator;

protected:
  std::string description;
  std::list<TestCallerInterface *> tests;

public:
  TestSuite(const std::string &desc);
  virtual ~TestSuite();

  void addTest(TestCallerInterface *test);

  const std::string &getDescription();

  iterator begin();
  iterator end();
};


class TestRunner
{
protected:
  std::ostream &stream;
  std::list<TestSuite *> suites;

public:
  TestRunner(std::ostream &stream);
  virtual ~TestRunner();

  void addSuite(TestSuite *suite);

  void operator() ();
};


#define UT_ASSERT(t) this->assertTrue(t, __FILE__, __LINE__)
#define UT_ASSERT_EQUAL(t, e) this->assertEqual(t, e, __FILE__, __LINE__)
#define UT_ASSERT_NEAR(t, e) this->assertNear(t, e, __FILE__, __LINE__)
#define UT_ASSERT_THROW(t, e) \
  try { t; throw UnitTest::TestFailure("No exception thrown!"); } catch (e &err) {}
}


#endif // UNITTEST_HH
