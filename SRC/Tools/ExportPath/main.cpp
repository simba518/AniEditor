#include <eigen3/Eigen/Dense>
#include <Objmesh.h>
#include <JsonFilePaser.h>
#include <AuxTools.h>
using namespace Eigen;
using namespace UTILITY;

int main(int argc, char *argv[]){

  const string d = "/home/simba/Workspace/AnimationEditor/Data/flower_box_casa/";
  const string ball_f = d+"model/ball.obj";
  Objmesh ball;
  bool succ = ball.load(ball_f);
  ERROR_LOG_COND("failed to load " << ball_f,succ);

  
  const string path_f = d+"model/circle_path.txt";
  JsonFilePaser json_f;
  succ = json_f.open(path_f);
  ERROR_LOG_COND("failed to load "<< path_f,succ);


  vector<int> con_frames;
  json_f.read("con_frames",con_frames);
  vector<double> con_pos;
  json_f.read("con_pos",con_pos);
  for (int i = 0; i < con_pos.size(); ++i){
    cout << con_pos[i] << ", ";
  }
  
  
  const Eigen::VectorXd vertices = ball.getVerts();
  Eigen::Vector3d cen = ball.getCenter();
  for (int i = 0; i < con_frames.size(); ++i){

    Vector3d p = -cen;
	p[0] += con_pos[i*3+0];
	p[1] += con_pos[i*3+1];
	p[2] += con_pos[i*3+2];
	
	VectorXd pos = vertices;
	for (int i = 0; i < pos.size()/3; i++){
	  pos[i*3+0] += p[0];
	  pos[i*3+1] += p[1];
	  pos[i*3+2] += p[2];
	}
	ball.setVerts(pos);
	string zeros = "";
	if (con_frames[i]<10){
	  zeros = "00";
	}else if (con_frames[i]<100){
	  zeros = "0";
	}

	const string save_to = d+"/tempt/ball"+zeros+TOSTR(con_frames[i])+".obj";
	cout << "save to : " << save_to << endl;
	succ = ball.write(save_to);
	ERROR_LOG_COND("failed to write: " << save_to, succ);
  }

  return 0;
}
