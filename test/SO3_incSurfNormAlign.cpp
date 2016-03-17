
#include <iostream>
#include <Eigen/Dense>
#include <manifold/SO3.h>
#include <manifold/S.h>

int main (int argc, char** argv) {

  uint32_t N = 10;
  uint32_t K = 2;
  
  double theta = -1.*M_PI/180.;
  Eigen::Matrix3d Rmu_;
  Rmu_ << 1, 0, 0,
         0, cos(theta), sin(theta),
         0, -sin(theta), cos(theta);
  SO3d Rmu(Rmu_);
  SO3d R = Rmu;
  double tau_R = 0.0000;

  std::vector<S3d> mus; 
  mus.push_back(S3d(Eigen::Vector3d(0.,cos(theta),-sin(theta))));
  mus.push_back(S3d(Eigen::Vector3d(0.,sin(theta),cos(theta))));
  std::vector<double> taus;
  taus.push_back(10.);
  taus.push_back(10.);

  theta = 45.*M_PI/180.;
  std::vector<S3d> ns; 
  std::vector<uint32_t> zs;
  for (uint32_t i=0; i<N/2; ++i) {
    ns.push_back(S3d(Eigen::Vector3d(0.,cos(theta),-sin(theta))));
    ns.push_back(S3d(Eigen::Vector3d(0.,sin(theta),cos(theta))));
    zs.push_back(0);
    zs.push_back(1);
  }
  
  double delta = 0.01;
  double f_prev = 1e99;
  double f = -tau_R*(Rmu.Inverse() + R).matrix().trace();
  for (uint32_t i=0; i<N; ++i)
    f -= taus[zs[i]]*mus[zs[i]].vector().transpose()*R.matrix()*ns[i].vector();
  std::cout << "f=" << f << std::endl;
  for (uint32_t it=0; it<100; ++it) {
    Eigen::Matrix3d J = -0.5*(Rmu.matrix() - (R+Rmu.Inverse()+R).matrix()); 
//    std::cout << J << std::endl;
    for (uint32_t i=0; i<N; ++i) {
      J -= 0.5*(taus[zs[i]]*(mus[zs[i]].vector()*ns[i].vector().transpose()
         - R.matrix()*(mus[zs[i]].vector()*ns[i].vector().transpose()*R.matrix())));
//      std::cout << mus[zs[i]].vector()*ns[i].vector().transpose()  << std::endl;
//      std::cout << mus[zs[i]].vector().transpose() << std::endl;
//      std::cout << ns[i].vector().transpose()  << std::endl;
    }
//    std::cout << J << std::endl;
    Eigen::Vector3d Jw = SO3d::vee(R.Inverse().matrix()*J);
    R += -delta*Jw;
//    std::cout << Jw.transpose() << std::endl;
    f_prev = f;
    f = -tau_R*(Rmu.Inverse() + R).matrix().trace();
    for (uint32_t i=0; i<N; ++i)
      f -= taus[zs[i]]*mus[zs[i]].vector().transpose()*R.matrix()*ns[i].vector();
//    if ((f_prev - f)/f < 1e-3) 
//      break;
    std::cout << "f=" << f << " df/f=" << (f_prev - f)/fabs(f)
//      << std::endl;
      << std::endl << R << std::endl;
  }
  std::cout << std::endl << Rmu << std::endl;
  std::cout << std::endl << R << std::endl;
  std::cout << acos(R.matrix()(1,1))*180/M_PI << std::endl;
}
