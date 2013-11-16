SET(APPNAME tempt2)

AUX_SOURCE_DIRECTORY(. test_src)

ADD_EXECUTABLE(${APPNAME} 
  ${test_src}
)

LINK_PROJECT(${APPNAME}

  editor
  warping
  basesim
  )

TARGET_LINK_LIBRARIES( ${APPNAME} 
  
  /usr/local/lib/libutility_ui.a
  /usr/local/lib/libutility_math.a
  /usr/local/lib/libutility_common.a
  /usr/local/lib/libutility_elastic.a
  ${Boost_LIBRARIES}
  ${GLEW_LIBRARY}
  ${OPENGL_LIBRARIES}

  ${ARPACK_LIBRARIES}
  ${SUPERLU2.0_LIBRARIES}
  ${CASADI_LIBRARIES}
  ${SUITESPARSE_LIBRARIES}
  ${IPOPT_LIBRARIES}
  /usr/lib/liblapack.a
  ${LAPACK_LIBRARIES}  
  libpthread.so
  )