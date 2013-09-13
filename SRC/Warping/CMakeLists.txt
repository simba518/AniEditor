SET(LIBRARY_NAME warping)

AUX_SOURCE_DIRECTORY(. warping_src)

ADD_LIBRARY(${LIBRARY_NAME}
  ${warping_src}
  )