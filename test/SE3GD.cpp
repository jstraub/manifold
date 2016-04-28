
#include <iostream>
#include <Eigen/Dense>
#include <manifold/SO3.h>
#include <manifold/gradientDescentSE3.h>

class GDSE3gmm : public GDSE3<double> {
 public:
  GDSE3gmm(double piA, const Eigen::Vector3d& muA, const
      Eigen::Matrix3d& covA, const Eigen::Vector3d& muB, const
      Eigen::Matrix3d& covB) 
    : piA_(piA), piB_(1.-piA), muA_(muA), muB_(muB), covA_(covA), covB_(covB)
  {};

  virtual void ComputeJacobian(const SE3d& theta, Eigen::Matrix<double,3,1>* J, double* f) {
    SE3d T = theta;
    Eigen::Matrix3d R = T.matrix().topLeftCorner(3,3);
    Eigen::Vector3d t = T.matrix().topRightCorner(3,1);
    Eigen::Vector3d m = R*muB_ - muA_;
    Eigen::Matrix3d SB = R*covB_*R.transpose();
    Eigen::Matrix3d S = covA_+SB;
    if (J) {


      Eigen::Matrix3d JMat = -0.5*((R.Inverse() + Rmu_).matrix() 
          - (Rmu_.Inverse() + R).matrix()); 
      *J = SO3d::vee(JMat);
//      std::cout << R << std::endl;
//      std::cout << Rmu_ << std::endl;
//      std::cout << JMat << std::endl;
//      std::cout << J->transpose() << std::endl;
    }
    if (f) {
      *f = -(Rmu_.Inverse() + R).matrix().trace();
    }
  };
 protected:
  SO3d Rmu_;
  double piA_;
  double piB_;
  Eigen::Vector3d muA_;
  Eigen::Vector3d muB_;
  Eigen::Matrix3d covA_;
  Eigen::Matrix3d covB_;
};

int main (int argc, char** argv) {
  
  SO3d R;
  std::cout << R << std::endl;

  double theta = 15.*M_PI/180.;
  Eigen::Matrix3d Rmu_;
  Rmu_ << 1, 0, 0,
         0, cos(theta), sin(theta),
         0, -sin(theta), cos(theta);
  SO3d Rmu(Rmu_);
  
  std::cout << Rmu << std::endl;
  std::cout << R+Rmu << std::endl;
  std::cout << R << std::endl;
  std::cout << Rmu+R << std::endl;

  std::cout << R-Rmu << std::endl;

  Eigen::Vector3d w = R-Rmu;
  std::cout << Rmu.Exp(w) << std::endl;

  std::cout << Rmu-R << std::endl;

  GDSO3vMF gd(Rmu);
  gd.Compute(R, 1e-8, 100);
  R = gd.GetMinimum();
  
//  double delta = 0.1;
//  double f_prev = 1e99;
//  double f = (Rmu.Inverse() + R).matrix().trace();
//  std::cout << "f=" << f << std::endl;
//  for (uint32_t it=0; it<100; ++it) {
//    Eigen::Matrix3d J = -0.5*((R.Inverse() + Rmu).matrix() - (Rmu.Inverse() + R).matrix()); 
//    Eigen::Vector3d Jw = SO3d::vee(J);
////    R = R.Exp(-delta*Jw);
//    R += -delta*Jw;
////    std::cout << Jw << std::endl;
//    f_prev = f;
//    f = (Rmu.Inverse() + R).matrix().trace();
////    if ((f_prev - f)/f < 1e-3) 
////      break;
//    std::cout << "f=" << f << " df/f=" << (f_prev - f)/f 
//      << std::endl;
////      << std::endl << R << std::endl;
//  }
  std::cout << std::endl << Rmu << std::endl;
  std::cout << std::endl << R << std::endl;
}
