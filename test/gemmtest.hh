/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef GEMMTEST_HH
#define GEMMTEST_HH

#include "unittest.hh"


class GEMMTest : public UnitTest::TestCase
{
public:
  void testRectRowMajor();
  void testRectTransposedRowMajor();
  void testSquareRowMajor();

  void testRectColMajor();
  void testRectTransposedColMajor();
  void testSquareColMajor();

public:
  static UnitTest::TestSuite *suite();
};

#endif // GEMMTEST_HH
