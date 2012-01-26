/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_BLAS_HH__
#define __LINALG_BLAS_HH__

#include "utils.hh"


/*
 * BLAS LEVEL 1 routoines:
 */
#include "dot.hh"
#include "nrm2.hh"

/*
 * BLAS LEVEL 2 routines.
 */
#include "gemv.hh"
#include "trmv.hh"
#include "getc2.hh"

/*
 * BLAS LEVEL 3 routines.
 */
#include "gemm.hh"
#include "trmm.hh"

#endif
