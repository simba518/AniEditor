SET(APPNAME EditorTest)

AUX_SOURCE_DIRECTORY(. test_src)

ADD_EXECUTABLE(${APPNAME} 
  ${test_src}
)

LINK_PROJECT(${APPNAME}

  editor
  warping
  basesim
  externallib
  basetestor
  )

TARGET_LINK_LIBRARIES( ${APPNAME} 

  /usr/local/lib/libutility_math.a
  /usr/local/lib/libutility_common.a
  ${Boost_LIBRARIES}
  ${GLM_LIBRARIES}
  ${GLEW_LIBRARY}
  ${OPENGL_LIBRARIES}

  ${CASADI_LIBRARIES}
  ${SUITESPARSE_LIBRARIES}
  ${IPOPT_LIBRARIES}
  /usr/lib/liblapack.a
  ${LAPACK_LIBRARIES}  
  libpthread.so
  )