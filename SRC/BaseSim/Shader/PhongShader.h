#ifndef _PHONGSHADER_H_
#define _PHONGSHADER_H_

#include <string>
#include <boost/shared_ptr.hpp>
#include "GLSLProgramObject.h"


namespace LSW_RENDERMESH{
  
  class PhongShader{
	
  public:
	static boost::shared_ptr<PhongShader> getInstance(){
	  
	  if(p_instance == NULL){
		p_instance = boost::shared_ptr<PhongShader>(new PhongShader());
	  }
	  return p_instance;
	}
	void initialize();
	void bind();
	void unbind();
	
  protected:
    PhongShader(){}
	
  private:
	CS::pGLSLProgramObject pSubFaceGLSL;
	static boost::shared_ptr<PhongShader> p_instance;
  };
  
  typedef boost::shared_ptr<PhongShader> pPhongShader;
  
}//end of namespace

#endif /*_PHONGSHADER_H_*/
