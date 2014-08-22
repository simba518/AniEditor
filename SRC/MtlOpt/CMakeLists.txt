SET(LIBRARY_NAME mtlopt)

AUX_SOURCE_DIRECTORY(. mtlopt_src)

INCLUDE_DIRECTORIES(
  /home/simba/mosek/7/tools/platform/linux64x86/h/
  )

ADD_LIBRARY(${LIBRARY_NAME}
  ${mtlopt_src}
  )