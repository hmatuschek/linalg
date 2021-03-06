/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_OPENMP_WRAPPER_HH__
#define __LINALG_OPENMP_WRAPPER_HH__

#include <cstdlib>
#include <omp.h>

#ifdef _OPENMP
#define LINALG_HAS_OPENMP
#endif

/**
 * Simple wrapper class to encapsulate OpenMP.
 */
class OpenMP
{
public:
  /**
   * Retruns the number of available threads.
   *
   * Returns 1, if OpenMP is disabled.
   */
  static inline size_t getMaxThreads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#endif
    return 1;
  }

  /**
   * Returns the thread-identifier or 0 if OpenMP is disabled.
   */
  static size_t getThreadNum() {
#ifdef _OPENMP
    return omp_get_thread_num();
#endif
    return 0;
  }
};

#endif
