#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <assertext.h>
#include <Log.h>
#include "FrameName.h"
#include "ObjFileSeqLoader.h"
using namespace std;
using namespace boost;
using namespace LSW_RENDERMESH;

bool ObjFileSeqLoader::load(const string &obj_directory, vector<pObjRenderMesh> &objs){

  TRACE_FUN();
  objs.clear();
  bool succ = true;

  const vector<string> unsort_files = getFileNames(obj_directory);
  const vector<string> files = sortFileNameByFrameNum(unsort_files);

  if(files.size() <= 0){
	WARN_LOG("no animation sequence will be loaded, as there is no .obj file in " << obj_directory);
	succ = false;
  }
  
  BOOST_FOREACH(const string &f, files){
	pObjRenderMesh obj_mesh = pObjRenderMesh(new ObjRenderMesh());
	DEBUG_LOG("begin to load " << f);
	if (!obj_mesh->read(f)){
	  succ = false;
	  ERROR_LOG("failed to load the .obj file from " << f);
	  break;
	}else{
	  objs.push_back(obj_mesh);
	}
  }

  return succ;
}

vector<string> ObjFileSeqLoader::getFileNames(const string &target_path){
  
  TRACE_FUN();
  vector<string> all_matching_files;

  if (!filesystem::exists(target_path)){
	ERROR_LOG("no such directory " << target_path);
	return all_matching_files;
  }

  const string obj_appendix = ".obj";
  filesystem::directory_iterator end_itr;
  for(filesystem::directory_iterator i( target_path ); i != end_itr; ++i ){
	
	if( !filesystem::is_regular_file( i->status() ) ){
	  continue;// skip if not a file 
	}

	const string file_path = i->path().filename().string();
	string::size_type pos = file_path.find_last_of(".");
	const string appendix = file_path.substr(pos,file_path.size());
	if( !(appendix == obj_appendix) ){
	  continue; // skip if no match 
	}
	
	// file matches, store it
	all_matching_files.push_back(target_path + "/" + file_path);
  }
  return all_matching_files;
}

vector<string> ObjFileSeqLoader::sortFileNameByFrameNum(const vector<string> &files){

  const int n = (int)files.size();
  vector<FrameName> tempt(n);
  vector<string> sorted_files(n);

  for (int i = 0; i < n; ++i){
    tempt[i].setFileName( files[i] );
  }
  
  sort(tempt.begin(),tempt.end());

  for (int i = 0; i < n; ++i){
    sorted_files[i] = tempt[i].toString();
  }

  return sorted_files;
}
