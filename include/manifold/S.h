#pragma once
#include <Eigen/Dense>
#include <manifold/manifold.h>

template<typename T, uint32_t D>
class S : M<T,D> {
 public:
  
  Eigen::Matrix<T,D,1> operator-(const S<T,D>& other);
 private:
};
