#include <Log.h>
#include "PhongShader.h"
using namespace LSW_RENDERMESH;

pPhongShader PhongShader::p_instance;

void PhongShader::initialize(){

  if (pSubFaceGLSL){
	// pSubFaceGLSL should be initialized only once.
	return ;
  }

  /// @todo it's better to read the path fome a file.
  const std::string home_dir = std::string(getenv("HOME"));
  const std::string SHADER_PATH = home_dir+"/Workspace/PathControl/SRC/BaseSim/Shader/shaders/";

  /// @note this function should be called after OpenGL initialized.
  const GLenum err = glewInit();
  if (err != GLEW_OK){
	ERROR_LOG("failed to initialize the glew." << glewGetErrorString(err));
  }else{
	if (!glewIsSupported("GL_VERSION_2_0 GL_VERSION_1_5 GL_ARB_multitexture GL_ARB_vertex_buffer_object")) {
	  ERROR_LOG("failed to initialize the glew, as certern API of glew is not supported");
	}else{
	  pSubFaceGLSL = CS::pGLSLProgramObject(new CS::GLSLProgramObject);
	  pSubFaceGLSL->attachVertexShader(SHADER_PATH+"subface_vertex.glsl");
	  pSubFaceGLSL->attachFragmentShader(SHADER_PATH+"subface_fragment.glsl");
	  pSubFaceGLSL->link();
	}
  }
}

void PhongShader::bind(){

  if (pSubFaceGLSL)
  	pSubFaceGLSL->bind();
}

void PhongShader::unbind(){

  if (pSubFaceGLSL)
  	pSubFaceGLSL->unbind();
}
