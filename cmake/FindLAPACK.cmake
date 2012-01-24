#
# Try to find LAPACK libraries
#

find_library(LAPACK_LIBRARY NAMES lapack)

find_library(BLAS_LIBRARY NAMES blas)

set(LAPACK_LIBRARIES ${LAPACK_LIBRARY} ${BLAS_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LAPACK  DEFAULT_MSG
                                  LAPACK_LIBRARY BLAS_LIBRARY)

mark_as_advanced(LAPACK_LIBRARY LAPACK_LIBRARIES)
