#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <manifold/manifold.h>

template<typename T, uint32_t D, Manifold<T,D> M>
class GD {
 public:
  GD() {};
  virtual ~GD() {};

  virtual void Compute(const M& theta0, T thr, uint32_t itMax);

  virtual void ComputeJacobian(Eigen::Matrix<T,D,1>* J, T* f) = 0;

  const M& GetState() {return theta_;}
 protected:
  T c_;
  T t_;
  M theta_;
  void LineSearch(Eigen::Matrix<T,D,1>* J, T* f);
};

GD::GD() : 
  c_(0.1), t_(0.3)
{}

void LineSearch(Eigen::Matrix<T,D,1>* J, T* f) {
  T delta = 1.;
  ComputeJacobian(J, f);
  T fNew = *f;
  M thetaNew = theta_;
  Eigen::Matrix<T,D,1> d = - (*J)/J->norm();
  T m = J->dot(d);
  while (f-fNew < -c_*m*delta && delta > 1e-16) {
    thetaNew += delta*d;
    ComputeJacobian(NULL, &fNew);
    std::cout << f-fNew << " <? " << -c_*m*delta 
      << " fNew=" << fNew << " delta=" << delta << std::endl;
    delta *= t_;
  }
  *J = delta*d;
  *f = fNew;
}

void Compute(const M& theta0, T thr, uint32_t itMax) {
  theta = theta0;
  M thetaPrev = theta0;
  Eigen::Matrix<T,D,1> J; 
  J.fill(0.);
  T fPrev = 1e12;
  T f = 1e10;
//  T delta = 1e-2;
  uint32_t it=0;
  while((fPrev-f)/fabs(f) > thr && it < itMax) {
    fPrev = f;
    LineSearch(&J, &f);
//    ComputeJacobian(&J, &f);
    thetaPrev = theta_;
    theta_ += - delta*J;
    std::cout << "@" << it << " f=" << f 
      << " df/f=" << (fPrev-f)/fabs(f) << std::endl;
    ++it;
  }
  if (f > fPrev) {
    theta_ = thetaPrev;
  }
}
