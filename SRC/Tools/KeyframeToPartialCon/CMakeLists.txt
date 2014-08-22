SET(APP_NAME key2par)

AUX_SOURCE_DIRECTORY(. key2par_src)

ADD_EXECUTABLE(${APP_NAME} 
  
  ${key2par_src}
  )

TARGET_LINK_LIBRARIES(${APP_NAME}
  
  /usr/local/lib/libutility_math.a
  /usr/local/lib/libutility_common.a
  ${Boost_LIBRARIES}
  )