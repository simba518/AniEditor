SET(LIBRARY_NAME basesim)

AUX_SOURCE_DIRECTORY(./ObjMesh objmesh_src)
AUX_SOURCE_DIRECTORY(./ConAndDrag conanddrag_src)
AUX_SOURCE_DIRECTORY(./ basesim_src)

ADD_LIBRARY(${LIBRARY_NAME}

  ${objmesh_src}
  ${conanddrag_src}
  ${basesim_src}
  ${shader_src}
  )