/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef GEMVTEST_HH
#define GEMVTEST_HH

#include "unittest.hh"


class GEMVTest : public UnitTest::TestCase
{
public:
  void testSquareRowMajor();
  void testSquareColMajor();
  void testSquareTransposedRowMajor();
  void testSquareTransposedColMajor();

  void testRectRowMajor();
  void testRectTransposedRowMajor();
  void testRectColMajor();
  void testRectTransposedColMajor();

public:
  static UnitTest::TestSuite *suite();
};

#endif // GEMVTEST_HH
