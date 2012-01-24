#include "unittest.hh"
#include <iostream>

#include "nrm2test.hh"
#include "dottest.hh"
#include "gemvtest.hh"
#include "gemmtest.hh"


using namespace UnitTest;


int main(void)
{
  // Construct test-runner
  TestRunner runner(std::cout);

  runner.addSuite(NRM2Test::suite());
  runner.addSuite(DOTTest::suite());
  runner.addSuite(GEMVTest::suite());
  runner.addSuite(GEMMTest::suite());

  runner();

  return 0;
};
