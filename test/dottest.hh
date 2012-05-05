/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef DOTTEST_HH
#define DOTTEST_HH

#include "unittest.hh"


class DOTTest : public UnitTest::TestCase
{
public:
  void testIncrFloat();
  void testIncrDouble();

  void testDenseFloat();
  void testDenseDouble();

  void testHugeIncr();
  void testHugeDense();


public:
  static UnitTest::TestSuite *suite();
};

#endif // DOTTEST_HH
