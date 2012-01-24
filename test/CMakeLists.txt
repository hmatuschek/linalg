#
# Build script for unit-tests:
#

SET(LINALG_TEST_SOURCES main.cc
    unittest.cc cputime.cc
    nrm2test.cc dottest.cc gemvtest.cc gemmtest.cc)
SET(LINALG_TEST_HEADERS
    unittest.hh cputime.hh
    nrm2test.hh dottest.hh gemvtest.hh gemmtest.hh)


ADD_EXECUTABLE(linalg-test ${LINALG_TEST_SOURCES})
TARGET_LINK_LIBRARIES(linalg-test ${LAPACK_LIBRARY} ${BLAS_LIBRARY} m)
