SET(APPNAME Cubature)

AUX_SOURCE_DIRECTORY(. src)
AUX_SOURCE_DIRECTORY(./cubacode src_2)

ADD_EXECUTABLE(${APPNAME} 
  ${src_2}
  ${src})

LINK_PROJECT(${APPNAME}

  editor
  warping
  basesim
  externallib
)

TARGET_LINK_LIBRARIES(${APPNAME}

  /usr/local/lib/libutility_ui.a
  /usr/local/lib/libutility_math.a
  /usr/local/lib/libutility_common.a
  ${Boost_LIBRARIES}
  ${JSONCPP_LIBRARIES}
  ${LOG4CPLUS_LIBRARIES}
  ${CASADI_LIBRARIES}
  ${SUITESPARSE_LIBRARIES}
  ${LAPACK_LIBRARIES}
)