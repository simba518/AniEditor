SET(APP_NAME ExportCubatures)

AUX_SOURCE_DIRECTORY(. ExportPath_src)

ADD_EXECUTABLE(${APP_NAME} 
  ${ExportPath_src}
  )

TARGET_LINK_LIBRARIES(${APP_NAME}
  
  /usr/local/lib/libutility_math.a
  /usr/local/lib/libutility_common.a
  /usr/local/lib/libutility_elastic.a
  ${Boost_LIBRARIES}
  ${ARPACK_LIBRARIES}
  ${SUPERLU2.0_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${SUITESPARSE_LIBRARIES}
  )