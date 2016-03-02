
template<typename T, int D>
S<T,D>::S() : x_(Eigen::Matrix<T,D,1>::Zero()) {
  x_(0) = 1.;
}

template<typename T, int D>
S<T,D>::S(const Eigen::Matrix<T,D,1>& x) : x_(x)  
{}

template<typename T, int D>
S<T,D>::S(const S<T,D>& other) : x_(other.x_)  
{}

template<typename T, int D>
std::ostream& operator<<(std::ostream& out, const S<T,D>& q) {
  out << q.vector().transpose();
  return out;
}

template<typename T, int D>
Eigen::Matrix<T,D,1> S<T,D>::operator-(const S<T,D>& other) {
  
}
