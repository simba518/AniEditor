SET(LIBRARY_NAME editor)
AUX_SOURCE_DIRECTORY(. editor_src)

ADD_LIBRARY(${LIBRARY_NAME}
  ${editor_src}
  )