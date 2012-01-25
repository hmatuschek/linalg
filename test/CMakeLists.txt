#
# Build script for unit-tests:
#

SET(LINALG_TEST_SOURCES main.cc
    unittest.cc cputime.cc
    matrixtest.cc trimatrixtest.cc nrm2test.cc dottest.cc gemvtest.cc gemmtest.cc trmmtest.cc)
SET(LINALG_TEST_HEADERS
    unittest.hh cputime.hh
    matrixtext.cc trimatrixtest.cc nrm2test.hh dottest.hh gemvtest.hh gemmtest.hh trmmtest.hh)


ADD_EXECUTABLE(linalg-test ${LINALG_TEST_SOURCES})
TARGET_LINK_LIBRARIES(linalg-test ${LAPACK_LIBRARY} ${BLAS_LIBRARY} m)
