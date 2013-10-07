#ifndef _MESHVTKIO_H_
#define _MESHVTKIO_H_

#include <eigen3/Eigen/Dense>
#include <assertext.h>
#include <Objmesh.h>
#include <TetMesh.h>
#include <IO.h>
using namespace Eigen;
using namespace UTILITY;

namespace LSW_SIM{

  /**
   * @class MeshVtkIO read and write mesh as vtk files.
   * 
   */
  class MeshVtkIO{
	
  public:
	static bool write(pTetMesh_const volMesh,const string filename){
	  if(volMesh)
		return write(*volMesh,filename);
	  return false;
	}
	static bool write(pTetMesh_const volMesh,const VectorXd &u,const string filename){
	  if(volMesh)
		return write(*volMesh,u,filename);
	  return false;
	}
	static bool write(const TetMesh &volMesh,const string filename){
	  const VectorXd u0 = VectorXd::Zero(volMesh.nodes().size()*3);
	  return write(volMesh,u0,filename);
	}
	static bool write(const TetMesh &volMesh,const VectorXd &u, const string filename){

	  assert_eq(u.size(),volMesh.nodes().size()*3);
	  vector<double> x;
	  volMesh.nodes(x);
	  VVec3d v(x.size()/3);
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

	static bool write(pObjmesh_const objMesh,const string filename){
	  if(objMesh)
		return write(*objMesh,filename);
	  return false;
	}
	static bool write(const Objmesh &objMesh,const string filename){

	  const VectorXd &v = objMesh.getVerts();
	  if(v.size() <= 0){
		return false;
	  }
	  VVec3d vertices(objMesh.getVertsNum());
	  for (int i = 0; i < objMesh.getVertsNum(); ++i){
		vertices[i][0] = v[i*3+0];
		vertices[i][1] = v[i*3+1];
		vertices[i][2] = v[i*3+2];
	  }

	  const VectorXi &tri = objMesh.getFaces();
	  VVec3i faces(tri.size()/3);
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
