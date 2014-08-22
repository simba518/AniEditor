#include <set>
#include <vector>
#include <string>
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <PartialConstraints.h>
#include <TetMesh.h>
#include <VTKWriter.h>
#include <MatrixIO.h>
#include <AuxTools.h>
using namespace std;
using namespace Eigen;
using namespace UTILITY;
using namespace EIGEN3EXT;

bool loadCubatures(const string cub_elements_str,const string cub_weights_str,
				   vector<std::pair<double, int> > &cub_weights_elements){
  
  vector<int> cub_elements;
  vector<double> cub_weights;
  bool succ = loadVec(cub_weights_str,cub_weights);
  succ &= loadVec(cub_elements_str,cub_elements,UTILITY::TEXT);
  if(succ){
	assert_eq(cub_elements.size(),cub_weights.size());
	for (size_t i =0; i < cub_weights.size(); ++i)
	  cub_weights_elements.push_back(make_pair(cub_weights[i],cub_elements[i]));
  }
  sort(cub_weights_elements.begin(),cub_weights_elements.end());
  cout<< "cubatures: " << cub_weights_elements.size() << endl;

  return succ;
}

int main(int argc, char *argv[]){
  
  TetMesh tetMesh;
  MatrixXd U;
  vector<std::pair<double, int> > cub_weights_elements;
  vector<int> keyframe_ids;
  int desired_points = -1;
  string save_partial_to;

  if (argc < 5){
	INFO_LOG("usage: tet_file keyframes_file cub_elements_file cub_weights_file keyframe_ids [num_desired_points]");
	return -1;
  }else{

	bool succ = tetMesh.load(argv[1]);
	ERROR_LOG_COND("failed to load tetrahedron file: "<<argv[1],succ);

	save_partial_to = string(argv[2])+"_parcon.txt";
	succ = EIGEN3EXT::load(argv[2],U);
	ERROR_LOG_COND("failed to load Uk: "<<argv[2],succ);

	succ = loadCubatures(argv[3],argv[4],cub_weights_elements);
	ERROR_LOG_COND("failed to load the cub_weights/cub_eles: "<<argv[3]<<"/"<<argv[4],succ);

	assert_ge(argc,5+U.cols());
	cout << "keyframe ids: ";
	for (int i = 5; i < 5+U.cols(); ++i){
	  const int f = atoi(argv[i]);
	  keyframe_ids.push_back(f);
	  cout << f << " ";
	}
	cout << endl;

	if (argc > 5+U.cols()+1)
	  desired_points = atoi(argv[5+U.cols()+1]);

	if (desired_points <=0)
	  desired_points = cub_weights_elements.size();
	cout << "desired points" <<  desired_points << endl;
  }

  // select constraint nodes 
  set<int> con_nodes_set;
  vector<int> con_nodes;
  const VVec4i &tets = tetMesh.tets();

  for (int i = 0; i < cub_weights_elements.size(); ++i){

	const int elem = cub_weights_elements[i].second;
	assert_in(elem,0,tets.size()-1);
	con_nodes_set.insert(tets[elem][0]);
	con_nodes_set.insert(tets[elem][1]);
	con_nodes_set.insert(tets[elem][2]);
	con_nodes_set.insert(tets[elem][3]);
	if (con_nodes_set.size() >= desired_points)
	  break;
  }
  
  BOOST_FOREACH(const int n, con_nodes_set){
	con_nodes.push_back(n);
  }
  cout << "positional constraints for each frame: " << con_nodes.size() << endl;

  // convert to partial constraints
  PartialConstraintsSet par_con;
  MatrixXd key_points(con_nodes.size()*3,U.cols());

  const VVec3d &nodes = tetMesh.nodes();
  assert_eq(keyframe_ids.size(),U.cols());
  for (size_t f = 0; f < U.cols(); ++f){

	const VectorXd& uf = U.col(f);
	Matrix<double,3,-1> pos_con(3,con_nodes.size());
	for (size_t i = 0; i < con_nodes.size(); ++i){

	  const int n = con_nodes[i];
	  assert_in(n*3,0,uf.size()-3);

	  pos_con(0,i) = uf(n*3+0);
	  pos_con(1,i) = uf(n*3+1);
	  pos_con(2,i) = uf(n*3+2);

	  key_points(i*3+0,f) = nodes[n][0]+pos_con(0,i);
	  key_points(i*3+1,f) = nodes[n][1]+pos_con(1,i);
	  key_points(i*3+2,f) = nodes[n][2]+pos_con(2,i);
	}
	par_con.addConNodes(con_nodes,keyframe_ids[f]);
	par_con.updatePc(pos_con,keyframe_ids[f]);
  }

  // save the partial constraints
  const bool succ = par_con.write(save_partial_to);
  ERROR_LOG_COND("failed to save the partial constraints "<<save_partial_to,succ);
  
  for (int f = 0; f < key_points.cols(); ++f){
    VTKWriter vtk;
	const VectorXd &p = key_points.col(f);
	vtk.addPoints(p);
	vtk.write(save_partial_to+TOSTR(f)+".vtk");
  }
  
  return 0;
}
