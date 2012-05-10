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
 * collection of functions directly interfacing the BLAS and LAPACK Fortran routines. This library
 * is intented to be easy usable and able to deal with NumPy arrays directly. This implies that
 * all methods provided by this library can handle matrices in any storage order. As some functions
 * of the LAPACK library are not capable to handle matrices in row-major storage, these functions
 * are reimplemented in C++ in an efficient way using SIMD instructions where it is possible.
 *
 * The most basic data type in this library is the @c Linalg::Array type, which represents a
 * multi-dimensional array of values. It is implemented as a class template to be able to deal with
 * any scalar value type. Note, that not all operations are implemented for all scalar types. There
 * are two linear algebra types derived from this array class, (a) the @c Linalg::Matrix class
 * representing a 2 dimensional array of scalars and (b) the @c Linalg::Vector type, representing
 * a column vector of scalars.
 *
 *
 */


/**
 * @defgroup matrix Matrix and vector types
 */



#ifndef __LINALG_HH__
#define __LINALG_HH__

#include "array.hh"
#include "vector.hh"
#include "matrix.hh"
#include "trimatrix.hh"
#include "symmatrix.hh"

#include "workspace.hh"
#include "exception.hh"

#include "blas/blas.hh"
#include "lapack/lapack.hh"


#endif
