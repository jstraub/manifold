
#include <iostream>
#include <Eigen/Dense>
#include <manifold/SO3.h>

int main (int argc, char** argv) {
  
  SO3d R;
  std::cout << R << std::endl;

  double theta = 45.*M_PI/180.;
  Eigen::Matrix3d R2_;
  R2_ << 1, 0, 0,
         0, cos(theta), sin(theta)
         0, -sin(theta), cos(theta);
  SO3d R2(R2_);
  
  std::cout << R2 << std::endl;
  std::cout << R+R2 << std::endl;
  std::cout << R2+R << std::endl;

  std::cout << R-R2 << std::endl;
  std::cout << R2-R << std::endl;

}
