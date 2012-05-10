/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_SSE_HH__
#define __LINALG_SSE_HH__

namespace Linalg {

/**
 * Template prototype of all SMID traits.
 */
template <class Scalar> class SIMDTraits;


/**
 * SMID traits for double vectors.
 */
template<>
class SIMDTraits<double>
{
public:
  /** Defines the elementary vector type for the scalar type. */
  typedef double vector __attribute__( (vector_size(16), aligned(8)) );

  /** Defines an union, that allows for a direct element access. */
  typedef union {
    /** The elements as vector type for operations. */
    vector   v;
    /** The elements as array for element access. */
    double   d[2];
  } uvector;

  /** Holds the number of elements in the vector/array. */
  const static size_t num_elements = 2;
};


/**
 * SMID traits for float vectors.
 */
template<>
class SIMDTraits<float>
{
public:
  typedef float vector __attribute__ ( (vector_size(16), aligned(4)) );

  typedef union {
    vector   v;
    float    d[4];
  } uvector;

  const static size_t num_elements = 4;
};


}

#endif // SSE_HH
