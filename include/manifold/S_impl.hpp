
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
  return other.Log(*this);  
}

template<typename T, int D>
S<T,D> S<T,D>::Exp(Eigen::Matrix<T,D,1>& x) {
  S<T,D> q;
  T theta = x.norm();
  if (fabs(theta) < 0.05)
  { // handle sin(0)/0
    q.vector() = p.vector()*cos(theta) + x*(1.-theta*theta/6.); // + O(x^4)
  }else{
    q.vector() = p.vector()*cos(theta) + x*(sin(theta)/theta);
  }
  return q;
}

template <typename T, int D>
T S<T,D>::invSincDot(T dot)
{
  // 2nd order taylor expansions for the limit cases obtained via mathematica
  if(static_cast<T>(MIN_DOT) < dot && dot < static_cast<T>(MAX_DOT))
    return acos(dot)/sqrt(1.-dot*dot);
  else if(dot <= static_cast<T>(MIN_DOT))
    return M_PI/(sqrt(2.)*sqrt(dot+1.)) -1. + M_PI*sqrt(dot+1.)/(4.*sqrt(2.))
      -(dot+1.)/3. + 3.*M_PI*(dot+1.)*sqrt(dot+1.)/(32.*sqrt(2.)) 
      - 2./15.*(dot+1.)*(dot+1.);
  else //if(dot >= static_cast<T>(MAX_DOT))
    return 1. - (dot-1)/3. + 2./5.*(dot-1.)*(dot-1.);
}

template<typename T, int D>
Eigen::Matrix<T,D,1> S<T,D>::Log(const S<T,D>& q) {
  T dot = max(static_cast<T>(-1.0),min(static_cast<T>(1.0),
        p.vector().dot(q.vector())));
  return (q.vector()-p.vector()*dot)*invSincDot(dot);
}

template<typename T, int D>
Eigen::Matrix<T,D-1,1> Intrinsic(const Eigen::Matrix<T,D,1>& x) {

  Eigen::Matrix<T,D,1> north;
  north.fill(0);
  north(D-1) = 1.;

  Matrix<T,D,D> northR = north_R_TpS2();
  Matrix<T,D,1> xhat = (northR * x);
  if(fabs(xhat(D_-1))>1e-5) 
  {
    // projection to zero last entry
    xhat -= xhat.dot(north)*north;
  }
  return xhat.topRows<D-1>();
}

template <typename T, int D>
Matrix<T,D,D> S<T,D>::north_R_TpS2() const
{
  Eigen::Matrix<T,D,1> north;
  north.fill(0);
  north(D-1) = 1.;
  return rotationFromAtoB<T>(p.vector(),north);
}

/* rotation from point A to B; percentage specifies how far the rotation will 
 * bring us towards B [0,1] */
template<typename T, int D>
Matrix<T,D,D> S<T,D>::rotationFromAtoB(const Matrix<T,D,1>& a, const
    Matrix<T,D,1>& b, T percentage=1.0)
{
  Matrix<T,D,D> bRa;
   
  T dot = b.transpose()*a;
//  ASSERT(fabs(dot) <=1.0, "a="<<a.transpose()<<" |.| "<<a.norm()
//      <<" b="<<b.transpose()<<" |.| "<<b.norm()
//      <<" -> "<<dot);
  dot = max(static_cast<T>(-1.0),min(static_cast<T>(1.0),dot));
//  cout << "dot="<<dot<<" | |"<<fabs(dot+1.)<<endl;
  if(fabs(dot -1.) < 1e-6)
  {
    // points are almost the same -> just put identity
    bRa =  Matrix<T,D,D>::Identity();
//    bRa(0,0) = cos(percentage*M_PI);
//    bRa(1,1) = cos(percentage*M_PI);
//    bRa(0,1) = -sin(percentage*M_PI);
//    bRa(1,0) = sin(percentage*M_PI);
  }else if(fabs(dot +1.) <1e-6) 
  {
    // direction does not matter since points are on opposing sides of sphere
    // -> pick one and rotate by percentage;
    bRa = -Matrix<T,D,D>::Identity();
    bRa(0,0) = cos(percentage*M_PI*0.5);
    bRa(1,1) = cos(percentage*M_PI*0.5);
    bRa(0,1) = -sin(percentage*M_PI*0.5);
    bRa(1,0) = sin(percentage*M_PI*0.5);
  }else{
    T alpha = acos(dot) * percentage;
//    cout << "alpha="<<alpha<<endl;

    Matrix<T,D,1> c;
    c = a - b*dot;
//    ASSERT(c.norm() >1e-5, "c="<<c.transpose()<<" |.| "<<c.norm());
    c /= c.norm();
    Matrix<T,D,D> A = b*c.transpose() - c*b.transpose();
    Matrix<T,D,D> temp = b*b.transpose() + c*c.transpose(); 
    T temp2 = cos(alpha)-1.; 
    bRa = Matrix<T,D,D>::Identity() + sin(alpha)*A + (temp2)*(temp);
  }
  return bRa;
}

