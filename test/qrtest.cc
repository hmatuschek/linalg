#include "qrtest.hh"
#include <qr.hh>

using namespace Linalg;


void
QRTest::setUp()
{
  // Allocate some space;
  tau = Vector<double>(3);

  static double a_data[9] = { 12, -51,   4,
                               6, 167, -68,
                              -4,  24, -41 };

  A = Matrix<double>::fromData(a_data, 3, 3, 3, 1, 0).copy();
}



void
QRTest::testRealNN()
{
  Matrix<double> B = A.copy();
  Vector<double> v(B.rows());
  geqrf(A, tau);

  std::cerr << A << std::endl;
  ormqr(A, tau, B, v, true, true);
  std::cerr << B << std::endl;

  for (size_t i=0; i<3; i++) {
    for (size_t j=i; j<3; j++) {
      UT_ASSERT_NEAR(B(i,j), A(i,j));
    }
  }
}



UnitTest::TestSuite *
QRTest::suite()
{
  UnitTest::TestSuite *s = new UnitTest::TestSuite("Tests for QR decomposition:");

  s->addTest(new UnitTest::TestCaller<QRTest>(
               "qr(double[n,n])", &QRTest::testRealNN));

  return s;
}
