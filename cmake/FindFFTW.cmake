# FFTW_INCLUDE_DIR = fftw3.h
# FFTW_LIBRARY = libfftw3.a
# FFTWF_LIBRARY = libfftw3f.a
# FFTW_FOUND = true if FFTW3 is found

find_path(FFTW_INCLUDE_DIR fftw3.h)
find_library(FFTW_LIBRARY fftw3)
find_library(FFTWF_LIBRARY fftw3f)

set(FFTW_FOUND FALSE)
if(FFTW_INCLUDE_DIR AND FFTW_LIBRARY AND FFTWF_LIBRARY)
    set(FFTW_FOUND TRUE)
endif()

mark_as_advanced(FFTW_INCLUDE_DIR FFTW_LIBRARY FFTWF_LIBRARY FFTW_FOUND)

MESSAGE("-- Found FFTW3:\n"
    "--   Include: ${FFTW_INCLUDE_DIR}\n"
    "--   Library: ${FFTW_LIBRARY}, ${FFTWF_LIBRARY}")