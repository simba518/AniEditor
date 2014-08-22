#include <iostream>
#include <Objmesh.h>
#include <TetMesh.h>
#include <MatrixIO.h>
#include <eigen3/Eigen/Dense>
using namespace std;
using namespace UTILITY;
using namespace Eigen;

double maxof3(const double v1, const double v2, const double v3){
  const double v = v1>v2 ? v1:v2;
  return v>v3 ? v:v3;
}

int main(int argc, char *argv[]){

  // load data
  const string d = "/home/simba/Workspace/AnimationEditor/Data/dino_scaled/";
  const string obj_file = d+"model/mesh.obj";
  Objmesh obj;
  obj.load(obj_file);

  TetMesh tet;
  const string tet_file = d+"model/mesh.abq";
  tet.load(tet_file);

  // get dimension
  const BBoxD bb = obj.getBBox();
  cout << "old bounding box: "<< endl;
  cout << bb.getWidth() << " "<< bb.getDeepth()<< " "<< bb.getHeight() << endl;
  const Vector3d &cen = bb.getCenter();
  const double max = maxof3(bb.getWidth(),bb.getDeepth(),bb.getHeight());
  assert_gt(max,0.0f);
  
  // scale obj
  Eigen::VectorXd v = obj.getVerts();
  for (int i = 0; i < v.size(); i+=3){
    v.segment<3>(i) -= cen;
  }
  v *= (1.0f/max);
  obj.setVerts(v);
  const string saveto = obj_file+"_scaled.obj";
  obj.write(saveto);
  const BBoxD bb2 = obj.getBBox();
  cout << "new bounding box: "<< endl;
  cout << bb2.getWidth() << " "<< bb2.getDeepth()<< " "<< bb2.getHeight() << endl;

  // scale tet
  VVec3d x = tet.nodes();
  for (int i = 0; i < x.size(); i++){
    x[i] -= cen;
	x[i] *= (1.0f/max);
  }
  tet.setRestPos(x);
  tet.write(tet_file+"_scaled.abq");
  tet.writeVTK(tet_file+"_scaled.vtk");

  // scale U
  MatrixXd U;
  EIGEN3EXT::load(d+"model/full_stvk_sim_continue.b",U);
  U *= (1.0f/max);
  const MatrixXd Usub = U.rightCols(U.cols()-10);
  EIGEN3EXT::write(d+"model/full_stvk_sim_continue_scaled.b",U);
  EIGEN3EXT::write(d+"model/full_stvk_sim_continue_begin10_scaled.b",Usub);
  tet.writeVTK(d+"tempt/mesh/scaled",U);

  return 0;
}
