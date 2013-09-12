SET(LIBRARY_NAME animationeditor)
AUX_SOURCE_DIRECTORY(. animationeditor_src)
AUX_SOURCE_DIRECTORY(./AD ad_src)
AUX_SOURCE_DIRECTORY(./Interpolators interpolators_src)
AUX_SOURCE_DIRECTORY(./Math math_src)

ADD_LIBRARY(${LIBRARY_NAME}  

  ${animationeditor_src}
  ${ad_src}
  ${interpolators_src}
  ${math_src}
  )