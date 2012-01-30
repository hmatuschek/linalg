#
# Collection of C++ templates to access Fortran functions.
#

SET(LINALG_BLAS_LEVEL1_HEADERS blas/dot.hh blas/nrm2.hh)
SET(LINALG_BLAS_LEVEL2_HEADERS blas/gemv.hh blas/getc2.hh blas/trmv.hh)
SET(LINALG_BLAS_LEVEL3_HEADERS blas/gemm.hh blas/trmm.hh blas/trsm.hh)
SET(LINALG_BLAS_HEADERS blas/blas.hh blas/utils.hh
    ${LINALG_BLAS_LEVEL1_HEADERS}
    ${LINALG_BLAS_LEVEL2_HEADERS}
    ${LINALG_BLAS_LEVEL3_HEADERS})

SET(LINALG_LAPACK_HEADERS lapack/lapack.hh
    lapack/trtrs.hh lapack/trtri.hh lapack/potrf.hh)

SET(LINALG_HEADERS
    linalg.hh memory.hh matrix.hh trimatrix.hh vector.hh exception.hh workspace.hh
    python.hh symmatrix.hh operators.hh)

SET(LINALG_SOURCES ${LINALG_HEADERS} ${LINALG_BLAS_HEADERS} ${LINALG_LAPACK_HEADERS})

ADD_CUSTOM_TARGET(linalg SOURCES ${LINALG_SOURCES})
