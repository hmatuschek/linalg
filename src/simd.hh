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

}

#endif // SSE_HH
