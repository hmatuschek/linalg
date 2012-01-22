
SET(LINALG_SOURCES
    array_base.cc blas.cc exception.cc lapack.cc workspace.cc)

SET(LINALG_HEADERS
    array_base.hh array_view.hh exception.hh lapack.hh linalg.hh workspace.hh)

ADD_LIBRARY(linalg SHARED ${LINALG_SOURCES})
TARGET_LINK_LIBRARIES(lapack blas m) 

