#include "choleskytest.hh"

using namespace Linalg;



void
CholeskyTest::setUp()
{
  static double a_data[9] = { 1,  2,  3,
                              2, 20, 26,
                              3, 26, 70 };

  static double a_chol_u_data[9] = { 1,2,3,
                                     0,4,5,
                                     0,0,6 };

  static double a_chol_l_data[9] = { 1,0,0,
                                     2,4,0,
                                     3,5,6 };

  A = Matrix<double>::fromData(a_data, 3, 3, 3, 1, 0).copy();
  Acu = Matrix<double>::fromData(a_chol_u_data, 3, 3, 3, 1, 0).copy();
  Acl = Matrix<double>::fromData(a_chol_l_data, 3, 3, 3, 1, 0).copy();


  static double b_data[18] = { 1,0, 2,-1,  3,-1,
                               2,1, 6, 0,  4,-2,
                               3,1, 4, 2, 29, 0 };

  static double b_chol_u_data[18] = {1,0, 2,-1,  3,-1,
                                     0,0, 1, 0, -3,-3,
                                     0,0, 0, 0,  1, 0 };

  static double b_chol_l_data[18] = { 1,0,  0,0, 0,0,
                                      2,1,  1,0, 0,0,
                                      3,1, -3,3, 1,0 };

  B = Matrix< std::complex<double> >::fromData((std::complex<double> *)b_data, 3, 3, 3, 1, 0).copy();
  Bcu = Matrix< std::complex<double> >::fromData((std::complex<double> *)b_chol_u_data, 3, 3, 3, 1, 0).copy();
  Bcl = Matrix< std::complex<double> >::fromData((std::complex<double> *)b_chol_l_data, 3, 3, 3, 1, 0).copy();
}


void
CholeskyTest::testRealLower()
{
  cholesky(A, false);
  // compare lower-tri matrix:
  for (size_t i=0; i<3; i++) {
    for (size_t j=0; j<=i; j++) {
      UT_ASSERT_EQUAL(A(i,j), Acl(i,j));
    }
  }
}

void
CholeskyTest::testRealUpper()
{
  cholesky(A, true);
  // compare upper-tri matrix:
  for (size_t i=0; i<3; i++) {
    for (size_t j=i; j<3; j++) {
      UT_ASSERT_EQUAL(A(i,j), Acu(i,j));
    }
  }
}


void
CholeskyTest::testCmplxLower()
{
  cholesky(B, false);
  // compare complex-tri matrix:
  for (size_t i=0; i<3; i++) {
    for (size_t j=0; j<=i; j++) {
      UT_ASSERT(B(i,j) == Bcl(i,j));
    }
  }
}

void
CholeskyTest::testCmplxUpper()
{
  cholesky(B, true);
  // compare complex-tri matrix:
  for (size_t i=0; i<3; i++) {
    for (size_t j=i; j<3; j++) {
      UT_ASSERT(B(i,j) == Bcu(i,j));
    }
  }
}



UnitTest::TestSuite *
CholeskyTest::suite()
{
  UnitTest::TestSuite *s = new UnitTest::TestSuite("Tests for cholesky()");

  s->addTest(new UnitTest::TestCaller<CholeskyTest>(
               "cholesky() (real-lower)",
               &CholeskyTest::testRealLower));

  s->addTest(new UnitTest::TestCaller<CholeskyTest>(
               "cholesky() (real-upper)",
               &CholeskyTest::testRealUpper));

  s->addTest(new UnitTest::TestCaller<CholeskyTest>(
               "cholesky() (complex-lower)",
               &CholeskyTest::testCmplxLower));

  s->addTest(new UnitTest::TestCaller<CholeskyTest>(
               "cholesky() (complex-upper)",
               &CholeskyTest::testCmplxUpper));

  return s;
}
