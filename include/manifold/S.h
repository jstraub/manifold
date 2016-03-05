#pragma once
#include <Eigen/Dense>
#include <manifold/manifold.h>

/// Class describing a data point on the sphere in D dimensions
template<typename T, int D>
class S : M<T,D> {
 public:
  S();
  S(const Eigen::Matrix<T,D,1>& x);
  S(const S<T,D>& other);
  ~S() {};

  Eigen::Matrix<T,D-1,1> Intrinsic(const Eigen::Matrix<T,D,1>& x);
//  S<T,D> Exp(Eigen::Matrix<T,D-1,1>& x) 
  S<T,D> Exp(Eigen::Matrix<T,D,1>& x);
  Eigen::Matrix<T,D,1> Log(const S<T,D>& q);
//  Eigen::Matrix<T,D-1,1> Log_north(const S<T,D>& q);

  const Eigen::Matrix<T,D,1>& vector() const {return p_;}
  Eigen::Matrix<T,D,1>& vector() {return p_;}
  
  Eigen::Matrix<T,D,1> operator-(const S<T,D>& other);

  Eigen::Matrix<T,D,D> north_R_TpS2() const;
  static Eigen::Matrix<T,D,D> rotationFromAtoB(const
      Eigen::Matrix<T,D,1>& a, const Eigen::Matrix<T,D,1>& b, T
      percentage=1.0);
  constexpr static double MIN_DOT=-0.98;
  constexpr static double MAX_DOT=0.98;
 private:
  Eigen::Matrix<T,D,1> p_;
  static T invSincDot(T dot);
};

typedef S<double,2> S2d;
typedef S<double,3> S3d;

template<typename T, int D>
std::ostream& operator<<(std::ostream& out, const S<T,D>& q);

#include <manifold/S_impl.hpp>
