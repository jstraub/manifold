#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <manifold/manifold.h>
#include <manifold/gradientDescent.h>

template <typename T>
class GDSO3 : public GD<T,3,SE3<T>> {
 public:
  GDSE3();
  ~GDSE3() {};
  virtual void ComputeJacobian(const SE3<T>& theta, Eigen::Matrix<T,6,1>* J, T* f) = 0;
 protected:
};

template <typename T>
GDSE3<T>::GDSE3() {
  this->c_=0.1;
  this->t_=0.3;
}
