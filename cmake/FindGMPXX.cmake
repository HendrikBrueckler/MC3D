# Try to find the GMPXX libraries
# GMPXX_FOUND - system has GMPXX lib
# GMPXX_INCLUDE_DIR - the GMPXX include directory
# GMPXX_LIBRARIES - Libraries needed to use GMPXX

# GMPXX needs GMP

find_package( GMP QUIET )

if(GMP_FOUND)

  if (GMPXX_INCLUDE_DIR AND GMPXX_LIBRARIES)
    # Already in cache, be silent
    set(GMPXX_FIND_QUIETLY TRUE)
  endif()

  find_path(GMPXX_INCLUDE_DIR NAMES gmpxx.h
            PATHS ${GMP_INCLUDE_DIR_SEARCH}
            DOC "The directory containing the GMPXX include files"
           )

  find_library(GMPXX_LIBRARY NAMES gmpxx
               PATHS ${GMP_LIBRARIES_DIR_SEARCH}
               DOC "Path to the GMPXX library"
               )

  # handle the QUIETLY and REQUIRED arguments and set GMPXX_FOUND to TRUE if
  # all listed variables are TRUE
  include(FindPackageHandleStandardArgs)
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(GMPXX
                          REQUIRED_VARS GMPXX_LIBRARY GMPXX_INCLUDE_DIR
                          VERSION_VAR GMPXX_VERSION_STRING)

  if(GMPXX_FOUND)
    set(GMPXX_LIBRARIES ${GMPXX_LIBRARY})
    set(GMPXX_INCLUDE_DIRS ${GMPXX_INCLUDE_DIR})
  endif()

  mark_as_advanced(GMPXX_INCLUDE_DIR GMPXX_LIBRARY)

endif()
