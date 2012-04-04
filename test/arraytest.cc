#include "arraytest.hh"

#include "array.hh"

using namespace Linalg;


void
ArrayTest::testAssignment()
{
  std::vector<size_t> dims(2); dims[0] = 3; dims[1] = 2;
  Array<double> A(dims), B;

  B = A;

  // A and B now share the same data:
  UT_ASSERT(B.ptr()==A.ptr());

  for(size_t i=0; i<3; i++) {
    for (size_t j=0; j<2; j++) {
      std::vector<size_t> idx(2); idx[0] = i; idx[1] = j;
      UT_ASSERT_EQUAL(B.at(idx), A.at(idx));
    }
  }
}


void
ArrayTest::testValueAssignment()
{
  std::vector<size_t> dims(2); dims[0] = 3, dims[1] = 2;
  Array<double> A(dims), B(dims);

  B.values() = A;

  // A and B do not share the same data:
  UT_ASSERT(B.ptr() != A.ptr());

  for(size_t i=0; i<3; i++) {
    for (size_t j=0; j<2; j++) {
      std::vector<size_t> idx(2); idx[0] = i; idx[1] = j;
      UT_ASSERT_EQUAL(B.at(idx), A.at(idx));
    }
  }

  // Test transposed assignment:
  dims[0] = 2; dims[1] = 3;
  B = Array<double>(dims);
  B.t().values() = A;

  for(size_t i=0; i<3; i++) {
    for (size_t j=0; j<2; j++) {
      std::vector<size_t> idx(2); idx[0] = i; idx[1] = j;
      UT_ASSERT_EQUAL(B.t().at(idx), A.at(idx));
    }
  }

}



UnitTest::TestSuite *
ArrayTest::suite()
{
  UnitTest::TestSuite *s = new UnitTest::TestSuite("Tests for Array class");

  s->addTest(new UnitTest::TestCaller<ArrayTest>(
               "double[m,n] assignment", &ArrayTest::testAssignment));

  s->addTest(new UnitTest::TestCaller<ArrayTest>(
               "double[m,n] value assignment", &ArrayTest::testValueAssignment));

  return s;
}
