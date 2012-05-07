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

/** Type of a two-element vector of doubles. */
typedef double __v2df __attribute__( (vector_size(16), aligned(8)) );

/** Type of a 4-element vector of flaots. */
typedef float  __v4sf __attribute__( (vector_size(16), aligned(4)) );

/** Type of a 2-element vector of doubles with element access. */
typedef union {
  __v2df   v;
  double   d[2];
} __vec2df;

/** Type of a 4-element vector of floats with element access. */
typedef union {
  __v4sf   v;
  float    d[4];
} __vec4sf;


/**
 * Template prototype of all SMID traits.
 */
template <class Scalar> class SMIDTraits;

/**
 * SMID traits for double vectors.
 */
template<>
class SMIDTraits<double>
{
  typedef double vector __attribute__ ( (vector_size(16), aligned(8)) );
  typedef union {
    vector   v;
    double   d[2];
  } uvector;

  const static size_t num_elements = 2;
};


/**
 * SMID traits for double vectors.
 */
template<>
class SMIDTraits<float>
{
  typedef float vector __attribute__ ( (vector_size(16), aligned(4)) );
  typedef union {
    vector   v;
    double   d[4];
  } uvector;

  const static size_t num_elements = 4;
};


}

#endif // SSE_HH
