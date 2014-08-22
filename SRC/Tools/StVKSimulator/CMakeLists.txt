SET(APPNAME StVKSimulator)

AUX_SOURCE_DIRECTORY(. tempt_src)

ADD_EXECUTABLE(${APPNAME} ${tempt_src})

TARGET_LINK_LIBRARIES( ${APPNAME} 
  /usr/local/lib/libutility_common.a
  /usr/local/lib/libutility_math.a
  /usr/local/lib/libutility_elastic.a
  ${Boost_LIBRARIES}
  ${SUITESPARSE_LIBRARIES}
  )
