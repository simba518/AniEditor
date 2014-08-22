SET(APP_NAME ExpandKeyToSeq)

AUX_SOURCE_DIRECTORY(. expandkeytoseq_src)

ADD_EXECUTABLE(${APP_NAME} 
  ${expandkeytoseq_src}
  )

TARGET_LINK_LIBRARIES(${APP_NAME}

  /usr/local/lib/libutility_common.a    
  ${Boost_LIBRARIES}
  )