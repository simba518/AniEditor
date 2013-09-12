SET(LIBRARY_NAME barbiclib)
AUX_SOURCE_DIRECTORY(volumetricMesh-v2.0 volumetricMesh_src)
AUX_SOURCE_DIRECTORY(./ barbiclib_src)

ADD_LIBRARY(${LIBRARY_NAME}
${volumetricMesh_src}
)