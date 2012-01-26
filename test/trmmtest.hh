/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef TRMMTEST_HH
#define TRMMTEST_HH

#include "unittest.hh"


class TRMMTest : public UnitTest::TestCase
{
public:
  void testUpperRowMajor();
  void testTransposedRowMajor();

  void testUpperColMajor();
  void testTransposedColMajor();

public:
  static UnitTest::TestSuite *suite();
};


#endif // TRMMTEST_HH
