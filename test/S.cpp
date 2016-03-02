#include <iostream>
#include <Eigen/Dense>
#include <manifold/S.h>

int main (int argc, char** argv) {
  
  S3d q;
  std::cout << q << std::endl;
  S3d mu;
  mu.vector() << 1./sqrt(2),1./sqrt(2), 0.;
  std::cout << mu << std::endl;

  double delta = 0.1;
  double f_prev = 1e99;
  double f = mu.vector().transpose()*q.vector();
  std::cout << "f=" << f << std::endl;
  for (uint32_t it=0; it<100; ++it) {
    Eigen::Vector3d J = -2.*(mu.vector() - q.vector()*q.vector().transpose()*mu.vector()); 
    q.vector() += -delta*J;
    q.vector() /= q.vector().norm();
//    std::cout << Jw << std::endl;
    f_prev = f;
    f = mu.vector().transpose()*q.vector();
//    if ((f_prev - f)/f < 1e-3) 
//      break;
    std::cout << "f=" << f << " df/f=" << (f_prev - f)/f 
      << std::endl;
  }
  std::cout << std::endl << mu << std::endl;
  std::cout << std::endl << q << std::endl;
}

