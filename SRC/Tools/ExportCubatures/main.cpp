#include <TetMesh.h>
#include <IOTools.h>
using namespace UTILITY;

int main(int argc, char *argv[]){

  // load tet mesh
  const string d = "/home/simba/Workspace/AnimationEditor/Data/flower_box/";
  const string tet_file = d + "/model/mesh.abq";
  TetMesh tetmesh;
  bool succ = tetmesh.load(tet_file);
  ERROR_LOG_COND("failed to load "<<tet_file, succ);
  

  // save surfaces as obj mesh
  const string save_tet_surface = d+"/tempt/volume_surface.obj";
  ofstream tet_obj(save_tet_surface);  
  ERROR_LOG_COND("failed to open: " << save_tet_surface, tet_obj.is_open());

  const VVec3d& nodes = tetmesh.nodes();
  const VVec3i& faces = tetmesh.surface();

  for (int i = 0; i < faces.size(); ++i){
	for (int j = 0; j < 3; ++j){
	  const int v = faces[i][j];
	  tet_obj <<"v "<<nodes[v][0]<<" "<<nodes[v][1]<<" "<<nodes[v][2]<<endl;
	}
	tet_obj << "f "<<i*3+1<<" "<<i*3+2<<" "<<i*3+3<<"\n";
  }

  // load cubatures
  const string cub_f = d+"/model/cubpoints.txt";
  vector<int> cub;
  succ = loadVec(cub_f,cub,UTILITY::TEXT);
  ERROR_LOG_COND("failed to load "<<cub_f, succ);

  // save cubatures
  const string save_cub_tet = d+"tempt/cub_tets.obj";
  ofstream cub_tet_boj(save_cub_tet);
  ERROR_LOG_COND("failed to open: " << save_cub_tet, cub_tet_boj.is_open());

  for (int i = 0; i < cub.size(); ++i){

    const tetrahedron& t = tetmesh.getTet(cub[i]);
	cub_tet_boj << "v "<<t._a(0,0)<<" "<<t._a(1,0)<<" "<<t._a(2,0)<<"\n";
	cub_tet_boj << "v "<<t._b(0,0)<<" "<<t._b(1,0)<<" "<<t._b(2,0)<<"\n";
	cub_tet_boj << "v "<<t._c(0,0)<<" "<<t._c(1,0)<<" "<<t._c(2,0)<<"\n";
	cub_tet_boj << "v "<<t._d(0,0)<<" "<<t._d(1,0)<<" "<<t._d(2,0)<<"\n";

	cub_tet_boj << "f "<<1+i*4+0<<" "<<1+i*4+2<<" "<<1+i*4+1<<"\n";
	cub_tet_boj << "f "<<1+i*4+0<<" "<<1+i*4+1<<" "<<1+i*4+3<<"\n";
	cub_tet_boj << "f "<<1+i*4+0<<" "<<1+i*4+3<<" "<<1+i*4+2<<"\n";
	cub_tet_boj << "f "<<1+i*4+1<<" "<<1+i*4+2<<" "<<1+i*4+3<<"\n";
  }
    
  return 0;
}
