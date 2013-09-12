#include <JsonFilePaser.h>
#include <Log.h>
#include <IEDSInterpolator.h>
#include <IMWInterpolator.h>
#include <RedRSInterpolator.h>
#include "InterpolatorFactory.h"
using namespace LSW_ANI_EDITOR;

pBaseInterpolator InterpolatorFactory::create(const string ini_file){
  
  pBaseInterpolator interpolator;
  UTILITY::JsonFilePaser inf;
  string interp_method;
  if ( inf.open(ini_file) && inf.read("interp_method", interp_method) ){
	if ( validMethod(interp_method) ){
	  interpolator = createMethod(interp_method);
	}
  }

  if ( interpolator != NULL && !interpolator->init (ini_file) ){
  	interpolator.reset();
  }
  return interpolator;
}

bool InterpolatorFactory::validMethod(const string method_name){

  const bool valid = ( string("IEDS") == method_name ||
					   string("IMW") == method_name||
					   string("REDRS") == method_name);
  if(!valid){
	ERROR_LOG("the method name is invalid!(in json node : interp_method )"<<
			  "the method name is " << method_name << ", which should be 'IEDS' or 'IMW' or 'REDRS'");
  }
  return valid;
}

pBaseInterpolator InterpolatorFactory::createMethod(const string in_method){
  
  pBaseInterpolator interpolator;
  if(string("IEDS") == in_method){
  	interpolator = pBaseInterpolator(new IEDSInterpolator());
  }else if(string("REDRS") == in_method){
  	interpolator = pBaseInterpolator(new RedRSInterpolator());
  }else{
  	interpolator = pIMWInterpolator(new IMWInterpolator());
  }
  return interpolator;
}
