#include "geqrftest.hh"
#include <lapack/geqrf.hh>
#include <lapack/ormqr.hh>

using namespace Linalg;


void
GEQRFTest::setUp()
{
  // Allocate some space;
  tau = Vector<double>(3);

  static double a_data[9] = { 12, -51,   4,
                               6, 167, -68,
                              -4,  24, -41 };

  A = Matrix<double>::fromData(a_data, 3, 3, 3, 1, 0).copy();
}



void
GEQRFTest::testRealNN()
{
  Matrix<double> B = A.copy();
  Vector<double> v(B.rows());
  Lapack::geqrf(A, tau);

  std::cerr << A << std::endl;
  Lapack::ormqr(A, tau, B, v, true, true);
  std::cerr << B << std::endl;

  for (size_t i=0; i<3; i++) {
    for (size_t j=i; j<3; j++) {
      UT_ASSERT_NEAR(B(i,j), A(i,j));
    }
  }
}



UnitTest::TestSuite *
GEQRFTest::suite()
{
  UnitTest::TestSuite *s = new UnitTest::TestSuite("Tests for geqrf :");

  s->addTest(new UnitTest::TestCaller<GEQRFTest>(
               "Lapack::geqrf(double[n,n])", &GEQRFTest::testRealNN));

  return s;
}
