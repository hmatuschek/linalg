#ifndef QRTEST_HH
#define QRTEST_HH

#include "unittest.hh"
#include <matrix.hh>
#include <vector.hh>


class QRTest : public UnitTest::TestCase
{
private:
  Linalg::Matrix<double> A;
  Linalg::Vector<double> tau;

public:
  virtual void setUp();

  void testRealNN();

public:
  static UnitTest::TestSuite *suite();
};

#endif // QRTEST_HH
