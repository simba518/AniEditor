SET(APP_NAME TetOpt)

ADD_EXECUTABLE(${APP_NAME} 
  ./main.cpp
  )

TARGET_LINK_LIBRARIES(${APP_NAME}
  /usr/local/lib/libutility_math.a
  /usr/local/lib/libutility_common.a  
  /usr/local/lib/libmesquite.a
  ${Boost_LIBRARIES}  
  )