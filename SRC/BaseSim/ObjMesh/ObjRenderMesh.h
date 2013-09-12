#ifndef OBJRENDERMESH_H
#define OBJRENDERMESH_H

#include <vector>
#include <string>
using std::vector;
using std::string;

#include <boost/shared_ptr.hpp>
#include <glm.h>

namespace LSW_RENDERMESH{
	
  class ObjRenderMesh{
	
  public:
	ObjRenderMesh();
	ObjRenderMesh(const string &file_name);
	~ObjRenderMesh();
	void setPhoneShading(bool use){
	  this->use_phone_shading = use;
	}
	bool usePhoneShading()const{
	  return use_phone_shading;
	}

	//file IO
	bool read(const string &file_name);
	bool save(const string &file_name,int mode = GLM_SMOOTH|GLM_TEXTURE|GLM_MATERIAL)const;

	//get information
	const GLMmodel *getModel()const{
	  return obj_model;
	}
	const string &filePath()const{
	  return file_path;
	}
	const float *Vertices()const;
	void Triangles(vector<int> &tri)const;
	int getVerticesNum()const;  // number of vertices 
	int getTriangleNum()const;  // number of triangles
	int getMaterialNum()const;  // number of materials for rendering
	int getTextureNum()const;   // number of textures
	
	void getDimension(float dimensions[3])const;
	void getCenter(double center_pos[3])const;

	/**
	 * mode is defined at <glm.h>, which could be:
	 * GLM_NONE     (0)             render with only vertices 
	 * GLM_FLAT     (1 << 0)        render with facet normals 
	 * GLM_SMOOTH   (1 << 1)        render with vertex normals 
	 * GLM_TEXTURE  (1 << 2)        render with texture coords 
	 * GLM_COLOR    (1 << 3)        render with colors 
	 * GLM_MATERIAL (1 << 4)        render with materials 
	 * GLM_2_SIDED  (1 << 5)        render two-sided polygons 
	 */
	void draw(const int mode)const;
	void draw()const;
	void draw(float color_r,float color_g,float color_b)const;
	int getRenderMode()const{
	  return default_render_mode;
	}

	//simple edit methods
	void scale(const float scale_f);
	void computeNormals();
	float unitize();
	void moveCenterTo(const double center_pos[3]);

	// set the vertices
	bool move(const vector<float> &u);
	bool setVertices(const float *vertices,int len);
	bool setVertices(const vector<float> &vertices);

	void reset();
	void clean();

  protected:
	void initDefaultRenderMode();

  protected:
	GLMmodel *obj_model;
	vector<float> ori_pos;
	int default_render_mode;
	string file_path;
	bool use_phone_shading; // default will use.
  };

  typedef boost::shared_ptr<ObjRenderMesh> pObjRenderMesh;
  typedef boost::shared_ptr<const ObjRenderMesh> pObjRenderMesh_const;
}

#endif
