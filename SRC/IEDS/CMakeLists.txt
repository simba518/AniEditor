SET(LIBRARY_NAME ieds)

AUX_SOURCE_DIRECTORY(. ieds_src)
AUX_SOURCE_DIRECTORY(./Warping warping_src)

ADD_LIBRARY(${LIBRARY_NAME}

  ${ieds_src}
  ${warping_src}
  )