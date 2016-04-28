#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <manifold/manifold.h>
#include <manifold/gradienDescent.h>

class GDSO3 : public GD<T,3,SO3<T>> {
 public:
  GDSO3();
  ~GDSO3() {};
  virtual void ComputeJacobian(Eigen::Matrix<T,D,1>* J, T* f) = 0;
 protected:
};

GDSO3::GDSO3() : 
  c_(0.1), t_(0.3)
{}
