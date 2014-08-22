#include <MatrixIO.h>
using namespace std;
using namespace Eigen;
using namespace EIGEN3EXT;

// expand one keyframe u to a sequence [u,u.....]
// ussage example: ExpandKeyToSeq ./keyframe_file.b 100
int main(int argc, char *argv[]){

  assert_ge_ext(argc,3,"ussage: ExpandKeyToSeq ./keyframe_file.b number_of_frames");
  const string key_file = argv[1];
  const int num_frames = atoi(argv[2]);
  VectorXd u;
  bool succ = load(key_file,u);
  assert_ext(succ, "failed to load: " << key_file);

  MatrixXd seq(u.size(), num_frames);
  for (int i = 0; i < seq.cols(); ++i)
    seq.col(i) = u;
  succ = write(key_file+".seq", seq);
  assert_ext(succ, "failed to write: " << key_file+".seq");

  return 0;
}
