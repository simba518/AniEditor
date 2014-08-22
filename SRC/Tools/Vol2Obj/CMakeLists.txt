SET(APP_NAME vol2obj)

AUX_SOURCE_DIRECTORY(. vol2obj_src)

ADD_EXECUTABLE(${APP_NAME} 
  ${vol2obj_src}
  )

TARGET_LINK_LIBRARIES(${APP_NAME}

  /usr/local/lib/libutility_common.a    
  ${Boost_LIBRARIES}
  )