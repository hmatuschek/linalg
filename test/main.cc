#include "unittest.hh"
#include <iostream>

#include "arraytest.hh"
#include "matrixtest.hh"
#include "trimatrixtest.hh"
#include "nrm2test.hh"
#include "dottest.hh"
#include "gemvtest.hh"
#include "trmvtest.hh"
#include "gemmtest.hh"
#include "trmmtest.hh"
#include "trsmtest.hh"

#include "trtrstest.hh"
#include "trtritest.hh"
#include "potrftest.hh"


using namespace UnitTest;


int main(void)
{
  // Construct test-runner
  TestRunner runner(std::cout);

  runner.addSuite(ArrayTest::suite());
  runner.addSuite(MatrixTest::suite());
  runner.addSuite(TriMatrixTest::suite());
  runner.addSuite(NRM2Test::suite());
  runner.addSuite(DOTTest::suite());
  runner.addSuite(GEMVTest::suite());
  runner.addSuite(TRMVTest::suite());
  runner.addSuite(GEMMTest::suite());
  runner.addSuite(TRMMTest::suite());
  runner.addSuite(TRSMTest::suite());

  runner.addSuite(TRTRSTest::suite());
  runner.addSuite(TRTRITest::suite());
  runner.addSuite(POTRFTest::suite());

  runner();

  return 0;
};
