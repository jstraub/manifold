#pragma once
#include <Eigen/Dense>
#include <manifold/manifold.h>

template<typename T, int D>
class S : M<T,D> {
 public:
  S();
  S(const Eigen::Matrix<T,D,1>& x);
  S(const S<T,D>& other);
  ~S() {};

  const Eigen::Matrix<T,D,1>& vector() const {return x_;}
  Eigen::Matrix<T,D,1>& vector() {return x_;}
  
  Eigen::Matrix<T,D,1> operator-(const S<T,D>& other);
 private:
  Eigen::Matrix<T,D,1> x_;
};

typedef S<double,2> S2d;
typedef S<double,3> S3d;

template<typename T, int D>
std::ostream& operator<<(std::ostream& out, const S<T,D>& q);

#include <manifold/S_impl.hpp>
