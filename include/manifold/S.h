#pragma once
#include <Eigen/Dense>
#include <manifold/manifold.h>
#include <manifold/SO3.h>

/// Class describing a data point on the sphere in D dimensions.
template<typename T, int D>
class S : M<T,D> {
 public:
  S();
  S(const Eigen::Matrix<T,D,1>& x);
  S(const S<T,D>& other);
  ~S() {};

  /// Compute the Riemannian Exp map around this point p of x in TpS.
  /// Yields another point on the sphere.
  S<T,D> Exp(const Eigen::Ref<const Eigen::Matrix<T,D,1>>& x) const;

  /// Compute the Riemannian Log map around this point p of q and
  /// yields a data point x in TpS. Note that the returned vector is
  /// represented in the ambient Euclidean space. See Intrinsic().
  Eigen::Matrix<T,D,1> Log(const S<T,D>& q) const;

  /// Compute the intrisic representation of a vector in TpS which is
  /// D-1 dimensional.
  Eigen::Matrix<T,D-1,1> Intrinsic(
      const Eigen::Ref<const Eigen::Matrix<T,D,1>>& x) const;

  /// Retraction that just orthogonally projects down to the sphere.
  /// A more efficient way of mapping TpS -> S than the Exp map.
  S<T,D> RetractOrtho(const Eigen::Ref<const Eigen::Matrix<T,D,1>>& x) const;

  /// Compute the dot product between two datapoints on the sphere.
  T dot(const S<T,D>& q) { return p_.dot(q.vector()); }

  /// Give access to the underlying vector.
  const Eigen::Matrix<T,D,1>& vector() const {return p_;}
  Eigen::Matrix<T,D,1>& vector() {return p_;}
  
  /// Compute the difference of two data points on S. The result will be
  /// in TotherS. This uses the Log map.
  Eigen::Matrix<T,D,1> operator-(const S<T,D>& other);

  /// Compute the rotation that rotates points in TpS to the north
  /// pole.
  Eigen::Matrix<T,D,D> north_R_TpS2() const;

  /// Compute the rotation that aligns a and b where a and b are on the
  /// sphere.
  static Eigen::Matrix<T,D,D> rotationFromAtoB(const
      Eigen::Matrix<T,D,1>& a, const Eigen::Matrix<T,D,1>& b, T
      percentage=1.0);

  constexpr static double MIN_DOT=-0.98;
  constexpr static double MAX_DOT=0.98;
 private:
  /// The data point.
  Eigen::Matrix<T,D,1> p_;
  
  /// Computes the inverse sinc(cos(dot)) in a stable way.
  static T invSincDot(T dot);
};

typedef S<double,2> S2d;
typedef S<double,3> S3d;

template<typename T, int D>
std::ostream& operator<<(std::ostream& out, const S<T,D>& q);

#include <manifold/S_impl.hpp>
