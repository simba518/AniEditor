#include <TetMesh.h>
using namespace UTILITY;

template <typename OS, typename INT>
void tet2vtk(OS &os,const double *node, size_t node_num,const INT *tet, size_t tet_num){

  os << "# vtk DataFile Version 2.0\nTET\nASCII\n\nDATASET UNSTRUCTURED_GRID\n";
  os << "POINTS " << node_num << " double\n";
  for(size_t i = 0; i < node_num; ++i)
	os << node[i*3+0] << " " << node[i*3+1] << " " << node[i*3+2] << "\n";

  os << "CELLS " << tet_num << " " << tet_num*5 << "\n";
  for(size_t i = 0; i < tet_num; ++i)
	os << 4 << "  "
	   << tet[i*4+0] << " " << tet[i*4+1] << " "
	   << tet[i*4+2] << " " << tet[i*4+3] << "\n";
  os << "CELL_TYPES " << tet_num << "\n";
  for(size_t i = 0; i < tet_num; ++i)
	os << 10 << "\n";
}

template <typename OS, typename Iterator, typename INT>
void cell_data(OS &os, Iterator first, INT size, const char *value_name, const char *table_name = "my_table"){
  os << "CELL_DATA " << size << "\n";
  vtk_data(os, first, size, value_name, table_name);
}

template <typename OS, typename Iterator, typename INT>
void vtk_data(OS &os, Iterator first, INT size, const char *value_name, const char *table_name = "my_table"){
  os << "SCALARS " << value_name << " double\nLOOKUP_TABLE " << table_name << "\n";
  for(size_t i = 0; i < size; ++i, ++first)
	os << *first << "\n";
}

void savemtl(const string modelname){
  
  const string d = "/home/simba/Workspace/AnimationEditor/Data/"+modelname+"/";
  const string tetfile = d+"/model/mesh.abq";
  const string mtlfile = d+"/model/mesh.elastic";

  TetMesh tetmesh;
  tetmesh.load(tetfile);
  tetmesh.loadElasticMtl(mtlfile);
  
  vector<double> nodes;
  tetmesh.nodes(nodes);
  vector<int> tets;
  tetmesh.tets(tets);

  vector<double> E;
  vector<double> rho;
  const ElasticMaterial<double> &mtl = tetmesh.material();
  for (int i = 0; i < tets.size()/4; ++i){
	rho.push_back(mtl._rho[i]);
	E.push_back(ElasticMaterial<double>::fromLameConstant(mtl._G[i],mtl._lambda[i])[0]);
  }
  
  const string f = d+"/tempt/mtl_E.vtk";
  ofstream file_E;
  file_E.open(f);


  tet2vtk(file_E, &(nodes[0]), (size_t)(nodes.size()/3), &(tets[0]), (size_t)(tets.size()/4));
  cell_data(file_E, E.begin(), E.size(), "Young's", "Young's");


  const string f2 = d+"/tempt/mtl_rho.vtk";
  ofstream file_rho;
  file_rho.open(f2);
  tet2vtk(file_rho, &(nodes[0]), (size_t)(nodes.size()/3), &(tets[0]), (size_t)(tets.size()/4));
  cell_data(file_rho, rho.begin(), rho.size(), "Density", "Density");
}

int main(int argc, char *argv[]){

  savemtl("dino_setmtl");
  savemtl("beam_setmtl");
  return 0;
}
