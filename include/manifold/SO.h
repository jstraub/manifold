#pragma once
#include <Eigen/Dense>
#include <manifold/manifold.h>

template<typename T>
class SO3 : M<T,3> {
 public:
  SO3();
  SO3(const Eigen::Matrix<T,3,3>& R);
  SO3(const SO3<T>& other);
  ~SO3() {};

  Eigen::Matrix<T,3,1> operator-(const SO3<T>& other);
  SO3<T>& operator+(const SO3<T>& other);

  SO3<T> Exp(const Eigen::Matrix<T,3,1>& w);
  Eigen::Matrix<T,3,1> Log(const SO3<T>& other);

  Eigen::Matrix<T,3,1> SO3<T,D>::vee()
  static Eigen::Matrix<T,3,3> SO3<T,D>::invVee(const Eigen::Matrix<T,3,1>& w);
 private:
  Eigen::Matrix<T,3,3> R_;

  Eigen::Matrix<T,D,1> vee(const Eigen::Matrix<T,3,3>& W);
  Eigen::Matrix<T,3,3> Exp(const Eigen::Matrix<T,3,1>& w);
  Eigen::Matrix<T,3,1> Log(const Eigen::Matrix<T,3,3>& R);
  
};

typedef SO3<double> SO3d;
typedef SO3<float> SO3f;

template<typename T>
std::ostream& operator<<(std::ostream& out, const SO3<T>& so3);
template<typename T>
SO3<T> operator+(SO3<T> lhs, const SO3<T>& rhs);

#include <manifold/SO3_impl.hpp>
