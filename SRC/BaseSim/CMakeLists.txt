SET(LIBRARY_NAME basesim)

AUX_SOURCE_DIRECTORY(./ObjMesh objmesh_src)
AUX_SOURCE_DIRECTORY(./ConAndDrag conanddrag_src)
AUX_SOURCE_DIRECTORY(./RedSimulator redsimulator_src)
AUX_SOURCE_DIRECTORY(./FullSimulator fullsimulator_src)
AUX_SOURCE_DIRECTORY(./Shader shader_src)
AUX_SOURCE_DIRECTORY(./ basesim_src)

ADD_LIBRARY(${LIBRARY_NAME}
${objmesh_src}
${conanddrag_src}
${basesim_src}
${redsimulator_src}
${fullsimulator_src}
${shader_src}
)