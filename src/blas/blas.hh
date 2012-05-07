/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

/**
 * @defgroup blas Interfaces to the BLAS library
 */


#ifndef __LINALG_BLAS_HH__
#define __LINALG_BLAS_HH__

#include "utils.hh"


/**
 * @defgroup blas1 BLAS level 1 routines
 * @ingroup blas
 */
#include "scal.hh"
#include "dot.hh"
#include "nrm2.hh"

/**
 * @defgroup blas2 BLAS level 2 routines
 * @ingroup blas
 */
#include "gemv.hh"
#include "trmv.hh"
#include "getc2.hh"

/**
 * @defgroup blas3 BLAS level 3 routines
 * @ingroup blas
 */
#include "gemm.hh"
#include "trmm.hh"
#include "trsm.hh"


/**
 * @defgroup blas_internal Internal used functions for BLAS routines.
 * @ingroup blas
 */

#endif
