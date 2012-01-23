#
# Collection of C++ templates to access Fortran functions.
#

SET(LINALG_BLAS_LEVEL1_HEADERS blas/dot.hh)
SET(LINALG_BLAS_LEVEL2_HEADERS blas/gemv.hh)
SET(LINALG_BLAS_LEVEL3_HEADERS blas/gemm.hh)
SET(LINALG_BLAS_HEADERS blas/blas.hh
    ${LINALG_BLAS_LEVEL1_HEADERS}
    ${LINALG_BLAS_LEVEL2_HEADERS}
    ${LINALG_BLAS_LEVEL3_HEADERS})

SET(LINALG_LAPACK_HEADERS lapack/lapack.hh lapack/getc2.hh)

SET(LINALG_HEADERS linalg.hh arraybase.hh matrix.hh vector.hh exception.hh workspace.hh)

SET(LINALG_SOURCES ${LINALG_HEADERS} ${LINALG_BLAS_HEADERS} ${LINALG_LAPACK_HEADERS})

ADD_CUSTOM_TARGET(linalg SOURCES ${LINALG_SOURCES})
