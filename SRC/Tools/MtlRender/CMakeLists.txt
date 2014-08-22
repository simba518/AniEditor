SET(APP_NAME MtlRender)

AUX_SOURCE_DIRECTORY(. KeyfApprox_src)

ADD_EXECUTABLE(${APP_NAME} 
  
  ${KeyfApprox_src}
  )

TARGET_LINK_LIBRARIES(${APP_NAME}
  
  /usr/local/lib/libutility_common.a
  ${Boost_LIBRARIES}
  )