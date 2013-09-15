#ifndef _MESHVTKIO_H_
#define _MESHVTKIO_H_

#include <boost/lexical_cast.hpp>
#include <assertext.h>
#include <eigen3/Eigen/Dense>
#include <VolObjMesh.h>
#include <ObjRenderMesh.h>
#include <IO.h>
using namespace Eigen;
using boost::lexical_cast;

namespace LSW_SIM{

  typedef Eigen::Vector3d Vec3;
  typedef Eigen::Matrix<sizeType, 3, 1> Vec3i;
  typedef std::vector<Vec3,Eigen::aligned_allocator<Vec3> >VectorV3;
  typedef vector<Vec3i,Eigen::aligned_allocator<Vec3i> >VectorV3i;
  
  /**
   * @class MeshVtkIO read and write mesh as vtk files.
   * 
   */
  class MeshVtkIO{
	
  public:
	static string int2str(const int number){
	  return lexical_cast<string>(number);
	}

	static bool write(pVolumetricMesh_const volMesh,const string filename){
	  if(volMesh)
		return write(*volMesh,filename);
	  return false;
	}
	static bool write(pVolumetricMesh_const volMesh,const VectorXd &u,const string filename){
	  if(volMesh)
		return write(*volMesh,u,filename);
	  return false;
	}
	static bool write(const VolumetricMesh &volMesh,const string filename){
	  const VectorXd u0 = VectorXd::Zero(volMesh.numVertices()*3);
	  return write(volMesh,u0,filename);
	}
	static bool write(const VolumetricMesh &volMesh,const VectorXd &u, const string filename){

	  assert_eq(u.size(),volMesh.numVertices()*3);
	  vector<double> x;
	  volMesh.vertices(x);
	  VectorV3 v(x.size()/3);
	  for (size_t i = 0; i < v.size(); ++i){
		v[i][0] = x[i*3+0] + u[i*3+0];
		v[i][1] = x[i*3+1] + u[i*3+1];
		v[i][2] = x[i*3+2] + u[i*3+2];
	  }
	  
	  COMMON::VTKWriter<double> writer("point",filename,true);
	  writer.appendPoints(v.begin(),v.end());
	  COMMON::VTKWriter<double>::IteratorIndex<COMMON::Vec3i> beg(0,0,1);
	  COMMON::VTKWriter<double>::IteratorIndex<COMMON::Vec3i> end(v.size(),0,1);
	  writer.appendCells(beg,end,COMMON::VTKWriter<double>::POINT);

	  return true;
	}

	static bool write(pObjRenderMesh_const objMesh,const string filename){
	  if(objMesh)
		return write(*objMesh,filename);
	  return false;
	}
	static bool write(const ObjRenderMesh &objMesh,const string filename){

	  const float *v = objMesh.Vertices();
	  if(!v){
		return false;
	  }
	  VectorV3 vertices(objMesh.getVerticesNum());
	  for (int i = 0; i < objMesh.getVerticesNum(); ++i){
		vertices[i][0] = v[i*3+0];
		vertices[i][1] = v[i*3+1];
		vertices[i][2] = v[i*3+2];
	  }

	  vector<int> tri;
	  objMesh.Triangles(tri);
	  VectorV3i faces(tri.size()/3);
	  for (int i = 0; i < faces.size(); ++i){
		faces[i][0] = tri[i*3+0];
		faces[i][1] = tri[i*3+1];
		faces[i][2] = tri[i*3+2];
	  }

	  COMMON::VTKWriter<double> writer("face",filename,true);
	  writer.appendPoints(vertices.begin(),vertices.end());
	  writer.appendCells(faces.begin(),faces.end(),COMMON::VTKWriter<double>::TRIANGLE);
	  return true;
	}
	
  };
  
}//end of namespace

#endif /*_MESHVTKIO_H_*/
