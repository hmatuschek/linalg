#
# Build script for unit-tests:
#

SET(BLAS1_TEST_SOURCES
    nrm2test.cc dottest.cc)
SET(BLAS1_TEST_HEADERS
    nrm2test.hh dottest.hh)

SET(BLAS2_TEST_SOURCES
    gemvtest.cc trmvtest.cc)
SET(BLAS2_TEST_HEADERS
    gemvtest.hh trmvtest.hh)

SET(BLAS3_TEST_SOURCES
    gemmtest.cc trmmtest.cc trsmtest.cc)
SET(BLAS3_TEST_HEADERS
    gemmtest.hh trmmtest.hh trsmtest.hh)

SET(LINALG_TEST_SOURCES main.cc
    unittest.cc cputime.cc matrixtest.cc trimatrixtest.cc
    ${BLAS1_TEST_SOURCES} ${BLAS2_TEST_SOURCES} ${BLAS3_TEST_SOURCES})
SET(LINALG_TEST_HEADERS
    unittest.hh cputime.hh matrixtext.hh trimatrixtest.hh
    ${BLAS1_TEST_HEADERS} ${BLAS2_TEST_HEADERS} ${BLAS3_TEST_HEADERS})


ADD_EXECUTABLE(linalg-test ${LINALG_TEST_SOURCES})
TARGET_LINK_LIBRARIES(linalg-test ${LAPACK_LIBRARY} ${BLAS_LIBRARY} m)
