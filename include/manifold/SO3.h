#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <manifold/manifold.h>

template<typename T>
class SO3 : M<T,3> {
 public:
  SO3();
  SO3(const Eigen::Matrix<T,3,3>& R);
  SO3(const SO3<T>& other);
  ~SO3() {};

  const Eigen::Matrix<T,3,3>& matrix() const { return R_;};
  Eigen::Matrix<T,3,3>& matrix() { return R_;};

  SO3<T> Exp (const Eigen::Matrix<T,3,1>& w)const ;
  Eigen::Matrix<T,3,1> Log (const SO3<T>& other)const ;

  Eigen::Matrix<T,3,1> vee() const;

  Eigen::Matrix<T,3,1> operator-(const SO3<T>& other);
  SO3<T>& operator+=(const SO3<T>& other);
  const SO3<T> operator+(const SO3<T>& other);

  static Eigen::Matrix<T,3,3> invVee(const Eigen::Matrix<T,3,1>& w);
  static Eigen::Matrix<T,3,1> vee(const Eigen::Matrix<T,3,3>& W);
  static Eigen::Matrix<T,3,3> skew(const Eigen::Matrix<T,3,3>& W);

 private:
  Eigen::Matrix<T,3,3> R_;

  static Eigen::Matrix<T,3,3> Exp_(const Eigen::Matrix<T,3,1>& w);
  static Eigen::Matrix<T,3,1> Log_(const Eigen::Matrix<T,3,3>& R);
  
};

typedef SO3<double> SO3d;
typedef SO3<float> SO3f;

template<typename T>
std::ostream& operator<<(std::ostream& out, const SO3<T>& so3);

#include <manifold/SO3_impl.hpp>
