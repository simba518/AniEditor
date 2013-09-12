# -Try to find SuperLU2.0
#
#
# The following are set after configuration is done: 
#  SuperLU_FOUND
#  SuperLU_LIBRARIES

SET(SUPERLU2.0_POSSIBLE_LIBPATHS

  /usr/local/lib
  /usr/lib
  /usr/lib64
  /usr/local/lib64
)

FIND_LIBRARY(SUPERLU2.0_LIBRARIES
  NAMES superlu2.0
  PATHS ${SUPERLU2.0_POSSIBLE_LIBPATHS}
)

#MESSAGE("DBG SUPERLU2.0_LIBBRARIES=${SUPERLU2.0_LIBRARIES})

IF(NOT SUPERLU2.0_LIBRARIES)
       MESSAGE(SEND_ERROR "SuperLU2.0 library not found.")
ENDIF(NOT SUPERLU2.0_LIBRARIES)

IF(SUPERLU2.0_LIBRARIES)
	MESSAGE(STATUS "Found SuperLU2.0: "${SUPERLU2.0_LIBRARIES})
	SET(SUPERLU2.0_FOUND TRUE)
ELSE(SUPERLU2.0_LIBRARIES)
	SET(SUPERLU2.0_FOUND FALSE)
ENDIF(SUPERLU2.0_LIBRARIES)

MARK_AS_ADVANCED(
  SUPERLU2.0_LIBRARIES
  SUPERLU2.0_FOUND
)