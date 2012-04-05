#include "trtritest.hh"

#include "matrix.hh"
#include "trimatrix.hh"
#include "lapack/trtri.hh"
#include "operators.hh"



using namespace Linalg;

void
TRTRITest::setUp()
{
  static double a_data[9] = { 1, 2, 3,
                       4, 5, 6,
                       7, 8, 9};

  static double bu_data[9] = { 1, -2./5, -1./15,
                        0,  1./5, -2./15,
                        0,     0,  1./9 };

  static double bl_data[9] =
  {     1.,      0,    0,
     -4./5,  1./5,    0,
    -1./15, -8./45, 1./9 };

  cA  = Matrix<double>::fromData( a_data, 3,3, 3,1, 0);
  cBu = Matrix<double>::fromData(bu_data, 3,3, 3,1, 0);
  cBl = Matrix<double>::fromData(bl_data, 3,3, 3,1, 0);

  fA  = Matrix<double>(3,3, false); fA.values() = cA;
  fBu = Matrix<double>(3,3, false); fBu.values() = cBu;
  fBl = Matrix<double>(3,3, false); fBl.values() = cBl;
}


void
TRTRITest::testUpperRowMajor()
{
  Matrix<double> tmp(cA.copy());
  Lapack::trtri(triu(tmp));

  for (size_t i=0; i<3; i++) {
    for (size_t j=i; j<3; j++) {
      UT_ASSERT_NEAR(tmp(i,j), cBu(i,j));
    }
  }
}


void
TRTRITest::testLowerRowMajor()
{
  Matrix<double> tmp(cA.copy());
  Lapack::trtri(tril(tmp));

  for (size_t i=0; i<3; i++) {
    for (size_t j=i; j<3; j++) {
      UT_ASSERT_NEAR(tmp(j,i), cBl(j,i));
    }
  }
}


void
TRTRITest::testUpperTransRowMajor()
{
  Matrix<double> tmp(cA.copy());
  Lapack::trtri(triu(tmp).t());

  for (size_t i=0; i<3; i++) {
    for (size_t j=i; j<3; j++) {
      UT_ASSERT_NEAR(tmp(i,j), cBu(i,j));
    }
  }
}


void
TRTRITest::testTransLowerRowMajor()
{
  Matrix<double> tmp(cA.copy());
  Lapack::trtri(tril(tmp.t()));

  for (size_t i=0; i<3; i++) {
    for (size_t j=i; j<3; j++) {
      UT_ASSERT_NEAR(tmp(i,j), fBu(i,j));
    }
  }
}



void
TRTRITest::testUpperColMajor()
{
  Matrix<double> tmp(fA.copy(false));
  Lapack::trtri(triu(tmp));

  for (size_t i=0; i<3; i++) {
    for (size_t j=i; j<3; j++) {
      UT_ASSERT_NEAR(tmp(i,j), fBu(i,j));
    }
  }
}



void
TRTRITest::testUpperTransColMajor()
{
  Matrix<double> A(4,4, false);

  A(0,0)= 1; A(0,1)= 2; A(0,2)= 3; A(0,3)= 4;
  A(1,0)= 5; A(1,1)= 6; A(1,2)= 7; A(1,3)= 8;
  A(2,0)= 9; A(2,1)=10; A(2,2)=11; A(2,3)=12;
  A(3,0)=13; A(3,1)=14; A(3,2)=15; A(3,3)=16;

  Matrix<double> tmp(A.copy());
  Lapack::trtri(triu(tmp).t());

  // Check upper triangular part for inverse:
  UT_ASSERT_NEAR(tmp(0,0), 1.);
  UT_ASSERT_NEAR(tmp(0,1), -1./3);
  UT_ASSERT_NEAR(tmp(0,2), -20./330);
  UT_ASSERT_NEAR(tmp(0,3), -375./9900);
  UT_ASSERT_NEAR(tmp(1,1), 5./30);
  UT_ASSERT_NEAR(tmp(1,2), -350./3300);
  UT_ASSERT_NEAR(tmp(1,3), -375./99000);
  UT_ASSERT_NEAR(tmp(2,2), 1./11);
  UT_ASSERT_NEAR(tmp(2,3), -675./9900);
  UT_ASSERT_NEAR(tmp(3,3), 625./10000);

  // Check lower part to be untouched:
  UT_ASSERT_EQUAL(tmp(1,0), A(1,0));
  UT_ASSERT_EQUAL(tmp(2,0), A(2,0));
  UT_ASSERT_EQUAL(tmp(2,1), A(2,1));
  UT_ASSERT_EQUAL(tmp(3,0), A(3,0));
  UT_ASSERT_EQUAL(tmp(3,1), A(3,1));
  UT_ASSERT_EQUAL(tmp(3,2), A(3,2));
}



void
TRTRITest::testTransLowerColMajor()
{
  Matrix<double> A(4,4, false);

  A(0,0)= 1; A(0,1)= 2; A(0,2)= 3; A(0,3)= 4;
  A(1,0)= 5; A(1,1)= 6; A(1,2)= 7; A(1,3)= 8;
  A(2,0)= 9; A(2,1)=10; A(2,2)=11; A(2,3)=12;
  A(3,0)=13; A(3,1)=14; A(3,2)=15; A(3,3)=16;

  Matrix<double> tmp(A.copy());
  Lapack::trtri(tril(tmp.t()));

  // Check upper triangular part for inverse:
  UT_ASSERT_NEAR(tmp(0,0), 1.);
  UT_ASSERT_NEAR(tmp(0,1), -1./3);
  UT_ASSERT_NEAR(tmp(0,2), -20./330);
  UT_ASSERT_NEAR(tmp(0,3), -375./9900);
  UT_ASSERT_NEAR(tmp(1,1), 5./30);
  UT_ASSERT_NEAR(tmp(1,2), -350./3300);
  UT_ASSERT_NEAR(tmp(1,3), -375./99000);
  UT_ASSERT_NEAR(tmp(2,2), 1./11);
  UT_ASSERT_NEAR(tmp(2,3), -675./9900);
  UT_ASSERT_NEAR(tmp(3,3), 625./10000);

  // Check lower part to be untouched:
  UT_ASSERT_EQUAL(tmp(1,0), A(1,0));
  UT_ASSERT_EQUAL(tmp(2,0), A(2,0));
  UT_ASSERT_EQUAL(tmp(2,1), A(2,1));
  UT_ASSERT_EQUAL(tmp(3,0), A(3,0));
  UT_ASSERT_EQUAL(tmp(3,1), A(3,1));
  UT_ASSERT_EQUAL(tmp(3,2), A(3,2));
}




UnitTest::TestSuite *
TRTRITest::suite()
{
  UnitTest::TestSuite *s = new UnitTest::TestSuite("Tests for Lapack::trtri()");

  s->addTest(new UnitTest::TestCaller<TRTRITest>(
               "Lapack::trtri(triu(double[m,m])) (row-major)", &TRTRITest::testUpperRowMajor));

  s->addTest(new UnitTest::TestCaller<TRTRITest>(
               "Lapack::trtri(tril(double[m,m])) (row-major)", &TRTRITest::testLowerRowMajor));

  s->addTest(new UnitTest::TestCaller<TRTRITest>(
               "Lapack::trtri(triu(double[m,m])::t()) (row-major)",
               &TRTRITest::testUpperTransRowMajor));

  s->addTest(new UnitTest::TestCaller<TRTRITest>(
               "Lapack::trtri(tril(double[m,m]::t())) (row-major)",
               &TRTRITest::testTransLowerRowMajor));

  s->addTest(new UnitTest::TestCaller<TRTRITest>(
               "Lapack::trtri(triu(double[m,m])) (col-major)", &TRTRITest::testUpperColMajor));

  s->addTest(new UnitTest::TestCaller<TRTRITest>(
               "Lapack::trtri(triu(double[m,m])::t()) (col-major)",
               &TRTRITest::testUpperTransColMajor));

  s->addTest(new UnitTest::TestCaller<TRTRITest>(
               "Lapack::trtri(tril(double[m,m]::t())) (col-major)",
               &TRTRITest::testTransLowerColMajor));

  return s;
}

