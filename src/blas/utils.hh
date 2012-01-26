/*
 * This file is part of the Linalg project, a C++ interface to BLAS and LAPACK.
 *
 * The source-code is licensed under the terms of the MIT license, read LICENSE for more details.
 *
 * (c) 2011, 2012 Hannes Matuschek <hmatuschek at gmail dot com>
 */

#ifndef __LINALG_BLAS_UTILS_HH__
#define __LINALG_BLAS_UTILS_HH__

/*
 * Some helper macros to simplify binding to Fortran code.
 */

#define BLAS_TO_COLUMN_MAJOR(A)     A.colMajor()
#define BLAS_ENSURE_COLUMN_MAJOR(A) (A = A.colMajor())

#define BLAS_IS_TRANSPOSED(A)       (A.isTransposed())
#define BLAS_IS_UPPER(T)            ((T.isUpper() && !T.isTransposed()) || (!T.isUpper() && T.isTransposed()))
#define BLAS_IS_LOWER(T)            ((!T.isUpper() && !T.isTransposed()) || (T.isUpper() && T.isTransposed()))
#define BLAS_HAS_UNIT_DIAG(T)       (T.hasUnitDiag())

#define BLAS_TRANSPOSED_FLAG(A)     BLAS_IS_TRANSPOSED(A) ? 'T' : 'N'
#define BLAS_UPLO_FLAG(T)           BLAS_IS_UPPER(T) ? 'U' : 'L'
#define BLAS_UNIT_DIAG_FLAG(T)      BLAS_HAS_UNIT_DIAG(T) ? 'U' : 'N'

#define BLAS_NUM_COLS(A)            (BLAS_IS_TRANSPOSED(A) ? A.rows() : A.cols())
#define BLAS_NUM_ROWS(A)            (BLAS_IS_TRANSPOSED(A) ? A.cols() : A.rows())
#define BLAS_LEADING_DIMENSION(A)   (A.stride())

#define BLAS_DIMENSION(x)           (x.dim())
#define BLAS_INCREMENT(x)           (x.stride())


#endif // __LINALG_BLAS_UTILS_HH__
