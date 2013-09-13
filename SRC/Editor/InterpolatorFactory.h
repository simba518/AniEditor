#ifndef _INTERPOLATORFACTORY_H_
#define _INTERPOLATORFACTORY_H_

#include <string>
#include <BaseInterpolator.h>
using namespace std;

namespace LSW_ANI_EDITOR{
  
  /**
   * @class InterpolatorFactory
   * 
   */
  class InterpolatorFactory{
	
  public:
	static pBaseInterpolator create(const string ini_file);
	static bool validMethod(const string method_name);
	static pBaseInterpolator createMethod(const string interp_method);

  };
  
}//end of namespace

#endif /*_INTERPOLATORFACTORY_H_*/
