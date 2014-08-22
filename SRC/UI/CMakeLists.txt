SET(APPNAME AniEditUI)

SET(LIBRARY_NAME anieditui) # compile as library for testing

SET(NONE_QT_CPP

  AniEditDM.cpp
  KeyframeDM.cpp
  AniEditDMConNodeDrag.cpp
)

#####################################QT#####################################
SET(UI_SOURCES

  VolApproxTri.cpp
  AniEditMainWin.cpp  
  AniEditDMRender.cpp
  AniEditDMConNodeDrag.cpp
)

SET(UI_HEADERS 

  ManipulateOP.h
  KeyframeDM.h
  KeyframeDMRender.h
  AniEditDM_UI.h
  AniEditMainWin.h
  AniEditDMRender.h
  AniEditDMConNodeDrag.h
  PreviewAniDMRender.h
  AniEditDMSelection.h
  KeyframeDMSelection.h
)

SET(UI_FORMS 

  aniedit.ui
  about.ui
)

QT4_WRAP_CPP(UI_HEADERS_MOC ${UI_HEADERS})
QT4_WRAP_UI(UI_FORMS_HEADERS ${UI_FORMS})

#####################################QT#####################################

ADD_EXECUTABLE(${APPNAME}

  # main_UITest.cpp  
  main.cpp
  ${NONE_QT_CPP}
  ${UI_SOURCES} 
  ${UI_HEADERS_MOC}
  ${UI_FORMS_HEADERS}
)

ADD_LIBRARY(${LIBRARY_NAME}

  ${NONE_QT_CPP}
  ${UI_SOURCES}
  ${UI_HEADERS_MOC}
  ${UI_FORMS_HEADERS}
)

LINK_PROJECT(${APPNAME}

  editor
  mtlopt
  warping
)

TARGET_LINK_LIBRARIES(${APPNAME}

  /usr/local/lib/libmosek64.so.7.0

  /usr/local/lib/libzjucadoptimizer.so
  /usr/local/lib/libhj-sparse-solver.so
  /usr/lib/mpich2/lib/libmpichcxx.so
 
  /usr/local/lib/libutility_ui.a
  /usr/local/lib/libutility_math.a
  /usr/local/lib/libutility_common.a
  /usr/local/lib/libutility_elastic.a
  /usr/local/lib/libalglib.a

  ${QGLViewer_LIBRARIES}
  ${QT_LIBRARIES}
  ${GLM_LIBRARIES}
  ${GLEW_LIBRARY}
  ${OPENGL_LIBRARIES}

  ${Boost_LIBRARIES}
  ${JSONCPP_LIBRARIES}
  ${LOG4CPLUS_LIBRARIES}

  ${ARPACK_LIBRARIES}
  ${SUPERLU2.0_LIBRARIES}

  ${SUITESPARSE_LIBRARIES}
  ${SUITESPARSE_SPQR_LIBRARY}
  ${ALGLIB_LIBRARIES}
  ${CASADI_LIBRARIES}
  ${IPOPT_LIBRARIES}
  ${LAPACK_LIBRARIES}
  # /usr/lib/liblapack.a
  # /usr/lib/libblas.so
  libpthread.so
)