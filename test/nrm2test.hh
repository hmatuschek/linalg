#ifndef NRM2TEST_HH
#define NRM2TEST_HH

#include "unittest.hh"


/**
 * Simple test-case for @c Linalg::BLAS::nrm2.
 */
class NRM2Test : public UnitTest::TestCase
{
public:
  void testVector();
  void testMatrixColumnRowMajor();
  void testMatrixColumnColMajor();
  void testMatrixRowRowMajor();
  void testMatrixRowColMajor();

public:
  static UnitTest::TestSuite *suite();
};

#endif // NRM2TEST_HH
