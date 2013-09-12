#ifndef _OBJFILESEQLOADER_H_
#define _OBJFILESEQLOADER_H_

#include <string>
#include <vector>
using namespace std;

#include <ObjMesh/ObjRenderMesh.h>
using namespace LSW_RENDERMESH;

namespace LSW_RENDERMESH{
  
  /**
   * @class ObjFileSeqLoader loader a set of obj meshes from a directory, which
   * are usually combine an animation sequence.
   * @see FrameName.h
   */
  class ObjFileSeqLoader{
	
  public:
	static bool load(const string &obj_directory, vector<pObjRenderMesh> &objs);
	
  protected: 
	static vector<string> getFileNames(const string &obj_directory);
	static vector<string> sortFileNameByFrameNum(const vector<string> &);

  };
  
}//end of namespace

#endif /*_OBJFILESEQLOADER_H_*/
