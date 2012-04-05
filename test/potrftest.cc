#include "potrftest.hh"
#include "lapack/potrf.hh"


using namespace Linalg;


void
POTRFTest::setUp()
{
  static double a_data[9] = { 1,  2,  3,
                              2, 20, 26,
                              3, 26, 70 };

  static double c_chol_u_data[9] = { 1,2,3,
                                     0,4,5,
                                     0,0,6 };

  cA = Matrix<double>::fromData(a_data, 3,3, 3,1, 0);
  fA = cA.copy(false);

  cACholU = Matrix<double>::fromData(c_chol_u_data, 3,3, 3,1, 0);
}


void
POTRFTest::testRowMajor()
{
  SymMatrix<double> tmp(cA.copy(true), true);
  Lapack::potrf(tmp);

  for(size_t i=0; i<3; i++) {
    for (size_t j=i; j<3; j++) {
      UT_ASSERT_NEAR(tmp(i,j), cACholU(i,j));
    }
  }
}


void
POTRFTest::testColMajor()
{
  SymMatrix<double> tmp(fA.copy(false), true);
  Lapack::potrf(tmp);

  for(size_t i=0; i<3; i++) {
    for (size_t j=i; j<3; j++) {
      UT_ASSERT_NEAR(tmp(i,j), cACholU(i,j));
    }
  }
}


UnitTest::TestSuite *
POTRFTest::suite()
{
  UnitTest::TestSuite *s = new UnitTest::TestSuite("Tests for Lapack::potrf()");

  s->addTest(new UnitTest::TestCaller<POTRFTest>(
               "Lapack::potrf(double[m,m]) (row-major)", &POTRFTest::testRowMajor));

  s->addTest(new UnitTest::TestCaller<POTRFTest>(
               "Lapack::potrf(double[m,m]) (col-major)", &POTRFTest::testColMajor));

  return s;
}
