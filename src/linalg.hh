/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

/**
 * @mainpage Linalg -- A C++ frontend to BLAS and LAPACK
 *
 * Linalg is a template library, that provides easy to use matrix and vector classes and a
 * collection of functions directly interfacing the BLAS and LAPACK Fortran routines.
 */


/**
 * @defgroup matrix Matrix and vector types
 */

/**
 * @defgroup operators Matrix and vector operations.
 * @ingroup matrix
 */


#ifndef __LINALG_HH__
#define __LINALG_HH__

#include "array.hh"
#include "vector.hh"
#include "matrix.hh"
#include "trimatrix.hh"

#include "workspace.hh"
#include "exception.hh"

#include "blas/blas.hh"
#include "lapack/lapack.hh"

#endif
