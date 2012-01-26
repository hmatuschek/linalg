/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef MATRIXTEST_HH
#define MATRIXTEST_HH

#include "unittest.hh"


class MatrixTest : public UnitTest::TestCase
{
public:
  void testRectRowToColMajor();
  void testSquareRowToColMajor();
  void testSubRowToColMajor();

  void testRectColToRowMajor();
  void testSquareColToRowMajor();
  void testSubColToRowMajor();

  void testFromData();
  void testValueRef();

  void testOwnershipTransfer();


public:
  static UnitTest::TestSuite *suite();
};

#endif // MATRIXTEST_HH
