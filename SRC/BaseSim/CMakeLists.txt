SET(LIBRARY_NAME basesim)

AUX_SOURCE_DIRECTORY(./ basesim_src)

ADD_LIBRARY(${LIBRARY_NAME}
  ${basesim_src}
  )