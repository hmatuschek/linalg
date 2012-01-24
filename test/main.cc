#include "unittest.hh"
#include <iostream>

#include "nrm2test.hh"
#include "dottest.hh"

using namespace UnitTest;


int main(void)
{
  // Construct test-runner
  TestRunner runner(std::cout);

  runner.addSuite(NRM2Test::suite());
  runner.addSuite(DOTTest::suite());

  runner();

  return 0;
};
