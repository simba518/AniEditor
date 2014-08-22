SET(APP_NAME unify)

AUX_SOURCE_DIRECTORY(. unify_src)

ADD_EXECUTABLE(${APP_NAME} 
  ${unify_src}
  )

TARGET_LINK_LIBRARIES(${APP_NAME}

  /usr/local/lib/libutility_common.a
  ${Boost_LIBRARIES}
  )