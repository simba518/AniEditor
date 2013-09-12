SET(APPNAME AniEditUI)

SET(LIBRARY_NAME anieditui) # compile as library for testing

SET(NONE_QT_CPP
  
  AniEditDM.cpp
  VolNodeGroupRender.cpp
  AniEditDMRender.cpp
  AniEditDMConNodeDrag.cpp
  AniEditDMVolVertexSel.cpp
  DragTrajectoryRecord.cpp
  DragGroupSelDraw.cpp
)

#####################################QT#####################################
SET(UI_SOURCES

  AniEditMainWin.cpp  
  VolObjMeshCtrl.cpp
)

SET(UI_HEADERS 

  AniEditDM_UI.h
  AniEditMainWin.h
  AniEditDMRenderCtrl.h
  AniEditDMSelCtrl.h
  AniEditDMDragCtrl.h
  PreviewAniDMRenderCtrl.h
  VolObjMeshCtrl.h
)

SET(UI_FORMS 

  aniedit.ui
  about.ui
)

QT4_WRAP_CPP(UI_HEADERS_MOC ${UI_HEADERS})
QT4_WRAP_UI(UI_FORMS_HEADERS ${UI_FORMS})

#####################################QT#####################################

ADD_EXECUTABLE(${APPNAME}

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

  animationeditor
  ieds
  basesim
  barbiclib
)

TARGET_LINK_LIBRARIES(${APPNAME}
 
  /usr/local/lib/libutility_ui.a
  ${QGLViewer_LIBRARIES}
  ${QT_LIBRARIES}
  ${GLM_LIBRARIES}
  ${GLEW_LIBRARY}
  ${OPENGL_LIBRARIES}

  ${Boost_LIBRARIES}
  ${JSONCPP_LIBRARIES}
  ${LOG4CPLUS_LIBRARIES}

  ${SUITESPARSE_LIBRARIES}
  ${SUITESPARSE_SPQR_LIBRARY}
  ${ALGLIB_LIBRARIES}
  ${CASADI_LIBRARIES}
  ${LAPACK_LIBRARIES}
)