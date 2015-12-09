#include <cmath>
#include <ostream>

/*
 * Written by: Jeffrey A. Fike
 * Stanford University, Department of Aeronautics and Astronautics
 *
 * Adapted to Eigen: Michael Tesch
 * Adapted to CppNumericalSolvers
 *
 * Copyright (c) 2006 Jeffrey A. Fike
 * Copyright (c) 2015 Michael Tesch
 * Copyright (c) 2015 Patrick Wieschollek
 */

template <typename eT>
class hyperdual {
  eT f0, f1, f2, f12;
 public:
  //creation operators and function to manually set values

  void wrt() {
    f1 = 1.;
    f2 = 0.;
    f12 = 0.;
  }
  void wrt(int i) {
    switch(i) {
      case 0:f0 = 1;break;
      case 1:f1 = 1;break;
      case 2:f2 = 1;break;
      case 3:f12 = 1;break;
    }
  }
  void nowrt() {
    f1 = 0.;
    f2 = 0.;
    f12 = 0.;
  }
  void nowrt(int i) {
    switch(i) {
      case 0:f0 = 0;break;
      case 1:f1 = 0;break;
      case 2:f2 = 0;break;
      case 3:f12 = 0;break;
    }
  }

  eT grad() {
    return f1;
  }
  eT hessian() {
    return f12;
  }

  hyperdual() {
    f0 = 0.0;
    f1 = 0.0;
    f2 = 0.0;
    f12 = 0.0;
  }

  hyperdual(eT x1, eT x2, eT x3, eT x4) {
    f0 = x1;
    f1 = x2;
    f2 = x3;
    f12 = x4;
  }

  hyperdual(eT x1) {
    f0 = x1;
    f1 = 0.0;
    f2 = 0.0;
    f12 = 0.0;
  }

  void setvalues(eT x1, eT x2, eT x3, eT x4) {
    f0 = x1;
    f1 = x2;
    f2 = x3;
    f12 = x4;
  }

  //examine values
  void view(void) {
    printf("%g  +  %g epsilon1  +  %g epsilon2  +  %g epsilon1 epsilon2\n",
           f0, f1, f2, f12);
  }

  eT real(void) const {
    return f0;
  }

  eT eps1(void) const {
    return f1;
  }

  eT eps2(void) const {
    return f2;
  }

  eT eps1eps2(void) const {
    return f12;
  }
  friend std::ostream & operator<<(std::ostream & output, const hyperdual & rhs) {
    // nolintnextline
    output << "{" << rhs.f0 << "," << rhs.f1 << "," << rhs.f2 << "," << rhs.f12 << "}";
    return output;
  }

  //basic manipulation
  hyperdual operator+() const {
    return *this;
  }
  hyperdual operator+(const hyperdual rhs) const {
    hyperdual temp;
    temp.f0 = f0 + rhs.f0;
    temp.f1 = f1 + rhs.f1;
    temp.f2 = f2 + rhs.f2;
    temp.f12 = f12 + rhs.f12;
    return temp;
  }

  friend hyperdual operator+(const eT lhs, const hyperdual rhs) {
    hyperdual temp;
    temp.f0 = lhs + rhs.f0;
    temp.f1 = rhs.f1;
    temp.f2 = rhs.f2;
    temp.f12 = rhs.f12;
    return temp;
  }
  hyperdual operator-() const {
    hyperdual temp;
    temp.f0 = -f0;
    temp.f1 = -f1;
    temp.f2 = -f2;
    temp.f12 = -f12;
    return temp;
  }
  hyperdual operator-(const hyperdual rhs) const {
    hyperdual temp;
    temp.f0 = f0 - rhs.f0;
    temp.f1 = f1 - rhs.f1;
    temp.f2 = f2 - rhs.f2;
    temp.f12 = f12 - rhs.f12;
    return temp;
  }
  friend hyperdual operator-(const eT lhs, const hyperdual rhs) {
    hyperdual temp;
    temp.f0 = lhs - rhs.f0;
    temp.f1 = -rhs.f1;
    temp.f2 = -rhs.f2;
    temp.f12 = -rhs.f12;
    return temp;
  }
  hyperdual operator*(const hyperdual rhs) const {
    hyperdual temp;
    temp.f0 = f0 * rhs.f0;
    temp.f1 = f0 * rhs.f1 + f1 * rhs.f0;
    temp.f2 = f0 * rhs.f2 + f2 * rhs.f0;
    temp.f12 = f0 * rhs.f12 + f1 * rhs.f2 + f2 * rhs.f1 + f12 * rhs.f0;
    return temp;
  }
  friend hyperdual operator*(const eT lhs, const hyperdual rhs) {
    hyperdual temp;
    temp.f0 = lhs * rhs.f0;
    temp.f1 = lhs * rhs.f1;
    temp.f2 = lhs * rhs.f2;
    temp.f12 = lhs * rhs.f12;
    return temp;
  }
  friend hyperdual operator/(const hyperdual lhs, const hyperdual rhs) {
    hyperdual temp, inv;
    inv = pow(rhs, -1);
    temp = lhs * inv;
    return temp;
  }
  friend hyperdual operator/(const eT lhs, const hyperdual rhs) {
    hyperdual temp, inv;
    inv = pow(rhs, -1);
    temp = lhs * inv;
    return temp;
  }
  friend hyperdual operator/(const hyperdual lhs, const eT rhs) {
    hyperdual temp;
    eT inv;
    inv = 1.0 / rhs;
    temp.f0 = inv * lhs.f0;
    temp.f1 = inv * lhs.f1;
    temp.f2 = inv * lhs.f2;
    temp.f12 = inv * lhs.f12;
    return temp;
  }

  friend eT ceil(const hyperdual lhs) {
    return ceil(lhs.f0);
  }

  hyperdual & operator+= (hyperdual rhs) {
    f0 += rhs.f0;
    f1 += rhs.f1;
    f2 += rhs.f2;
    f12 += rhs.f12;
    return *this;
  }
  hyperdual & operator-= (hyperdual rhs) {
    f0 -= rhs.f0;
    f1 -= rhs.f1;
    f2 -= rhs.f2;
    f12 -= rhs.f12;
    return *this;
  }
  hyperdual & operator*= (hyperdual rhs) {
    eT tf0, tf1, tf2, tf12;
    tf0 = f0;
    tf1 = f1;
    tf2 = f2;
    tf12 = f12;
    f0 = tf0 * rhs.f0;
    f1 = tf0 * rhs.f1 + tf1 * rhs.f0;
    f2 = tf0 * rhs.f2 + tf2 * rhs.f0;
    f12 = tf0 * rhs.f12 + tf1 * rhs.f2 + tf2 * rhs.f1 + tf12 * rhs.f0;
    return *this;
  }
  hyperdual & operator*= (eT rhs) {
    f0 *= rhs;
    f1 *= rhs;
    f2 *= rhs;
    f12 *= rhs;
    return *this;
  }

  hyperdual & operator/= (eT rhs) {
    f0 /= rhs;
    f1 /= rhs;
    f2 /= rhs;
    f12 /= rhs;
    return *this;
  }

  hyperdual & operator/= (hyperdual rhs) {
    *this = *this / rhs;
    return *this;
  }

  // operator eT() { return f0; }

  //math.h functions
  friend hyperdual pow(hyperdual x, eT a) {
    hyperdual temp;
    eT deriv, xval, tol;
    xval = x.f0;
    tol = 1e-15;
    if (fabs(xval) < tol) {
      if (xval >= 0)
        xval = tol;
      if (xval < 0)
        xval = -tol;
    }
    deriv = a * pow(xval, (a - 1));
    //temp.f0 = pow(xval,a);
    temp.f0 = pow(x.f0, a);       //Use actual x value, only use tol for derivs
    temp.f1 = x.f1 * deriv;
    temp.f2 = x.f2 * deriv;
    temp.f12 = x.f12 * deriv + a * (a - 1) * x.f1 * x.f2 * pow(xval, (a - 2));

    return temp;
  }
  friend hyperdual pow(hyperdual x, hyperdual a) {
    return exp(a * log(x));
  }
  friend hyperdual exp(hyperdual x) {
    hyperdual temp;
    eT deriv;
    deriv = exp(x.f0);
    temp.f0 = deriv;
    temp.f1 = deriv * x.f1;
    temp.f2 = deriv * x.f2;
    temp.f12 = deriv * (x.f12 + x.f1 * x.f2);
    return temp;
  }
  friend hyperdual log(hyperdual x) {
    hyperdual temp;
    eT deriv1, deriv2;
    deriv1 = x.f1 / x.f0;
    deriv2 = x.f2 / x.f0;
    temp.f0 = log(x.f0);
    temp.f1 = deriv1;
    temp.f2 = deriv2;
    temp.f12 = x.f12 / x.f0 - (deriv1 * deriv2);
    return temp;
  }
  friend hyperdual sin(hyperdual x) {
    hyperdual temp;
    eT funval, deriv;
    funval = sin(x.f0);
    deriv = cos(x.f0);
    temp.f0 = funval;
    temp.f1 = deriv * x.f1;
    temp.f2 = deriv * x.f2;
    temp.f12 = deriv * x.f12 - funval * x.f1 * x.f2;
    return temp;
  }
  friend hyperdual cos(hyperdual x) {
    hyperdual temp;
    eT funval, deriv;
    funval = cos(x.f0);
    deriv = -sin(x.f0);
    temp.f0 = funval;
    temp.f1 = deriv * x.f1;
    temp.f2 = deriv * x.f2;
    temp.f12 = deriv * x.f12 - funval * x.f1 * x.f2;
    return temp;
  }
  friend hyperdual tan(hyperdual x) {
    hyperdual temp;
    eT funval, deriv;
    funval = tan(x.f0);
    deriv = funval * funval + 1.0;
    temp.f0 = funval;
    temp.f1 = deriv * x.f1;
    temp.f2 = deriv * x.f2;
    temp.f12 = deriv * x.f12 + x.f1 * x.f2 * (2 * funval * deriv);
    return temp;
  }
  friend hyperdual asin(hyperdual x) {
    hyperdual temp;
    eT funval, deriv1, deriv;
    funval = asin(x.f0);
    deriv1 = 1.0 - x.f0 * x.f0;
    deriv = 1.0 / sqrt(deriv1);
    temp.f0 = funval;
    temp.f1 = deriv * x.f1;
    temp.f2 = deriv * x.f2;
    temp.f12 = deriv * x.f12 + x.f1 * x.f2 * (x.f0 * pow(deriv1, -1.5));
    return temp;
  }
  friend hyperdual acos(hyperdual x) {
    hyperdual temp;
    eT funval, deriv1, deriv;
    funval = acos(x.f0);
    deriv1 = 1.0 - x.f0 * x.f0;
    deriv = -1.0 / sqrt(deriv1);
    temp.f0 = funval;
    temp.f1 = deriv * x.f1;
    temp.f2 = deriv * x.f2;
    temp.f12 = deriv * x.f12 + x.f1 * x.f2 * (-x.f0 * pow(deriv1, -1.5));
    return temp;
  }
  friend hyperdual atan(hyperdual x) {
    hyperdual temp;
    eT funval, deriv1, deriv;
    funval = atan(x.f0);
    deriv1 = 1.0 + x.f0 * x.f0;
    deriv = 1.0 / deriv1;
    temp.f0 = funval;
    temp.f1 = deriv * x.f1;
    temp.f2 = deriv * x.f2;
    temp.f12 = deriv * x.f12 + x.f1 * x.f2 * (-2 * x.f0 / (deriv1 * deriv1));
    return temp;
  }
  friend hyperdual sqrt(hyperdual x) {
    return pow(x, 0.5);
  }
  friend hyperdual fabs(hyperdual x) {
    hyperdual temp;
    if (x < 0.0)
      temp = -x;
    else
      temp = x;
    return temp;
  }
  friend hyperdual max(hyperdual x1, hyperdual x2) {
    hyperdual temp;
    if (x1 > x2)
      temp = x1;
    else
      temp = x2;
    return temp;
  }
  friend hyperdual max(hyperdual x1, eT x2) {
    hyperdual temp;
    if (x1 > x2)
      temp = x1;
    else
      temp = x2;
    return temp;
  }
  friend hyperdual max(eT x1, hyperdual x2) {
    hyperdual temp;
    if (x1 > x2)
      temp = x1;
    else
      temp = x2;
    return temp;
  }
  friend hyperdual min(hyperdual x1, hyperdual x2) {
    hyperdual temp;
    if (x1 < x2)
      temp = x1;
    else
      temp = x2;
    return temp;
  }
  friend hyperdual min(hyperdual x1, eT x2) {
    hyperdual temp;
    if (x1 < x2)
      temp = x1;
    else
      temp = x2;
    return temp;
  }
  friend hyperdual min(eT x1, hyperdual x2) {
    hyperdual temp;
    if (x1 < x2)
      temp = x1;
    else
      temp = x2;
    return temp;
  }

  //comparisons
  friend bool operator>(hyperdual lhs, hyperdual rhs) {
    return (lhs.f0 > rhs.f0);
  }

  friend bool operator>(eT lhs, hyperdual rhs) {
    return (lhs > rhs.f0);
  }

  friend bool operator>(hyperdual lhs, eT rhs) {
    return (lhs.f0 > rhs);
  }

  friend bool operator>= (hyperdual lhs, hyperdual rhs) {
    return (lhs.f0 >= rhs.f0);
  }

  friend bool operator>= (eT lhs, hyperdual rhs) {
    return (lhs >= rhs.f0);
  }

  friend bool operator>= (hyperdual lhs, eT rhs) {
    return (lhs.f0 >= rhs);
  }

  friend bool operator<(hyperdual lhs, hyperdual rhs) {
    return (lhs.f0 < rhs.f0);
  }

  friend bool operator<(eT lhs, hyperdual rhs) {
    return (lhs < rhs.f0);
  }

  friend bool operator<(hyperdual lhs, eT rhs) {
    return (lhs.f0 < rhs);
  }

  friend bool operator<= (hyperdual lhs, hyperdual rhs) {
    return (lhs.f0 <= rhs.f0);
  }

  friend bool operator<= (eT lhs, hyperdual rhs) {
    return (lhs <= rhs.f0);
  }

  friend bool operator<= (hyperdual lhs, eT rhs) {
    return (lhs.f0 <= rhs);
  }

  friend bool operator== (hyperdual lhs, hyperdual rhs) {
    return (lhs.f0 == rhs.f0);
  }

  friend bool operator== (eT lhs, hyperdual rhs) {
    return (lhs == rhs.f0);
  }

  friend bool operator== (hyperdual lhs, eT rhs) {
    return (lhs.f0 == rhs);
  }

  friend bool operator!= (hyperdual lhs, hyperdual rhs) {
    return (lhs.f0 != rhs.f0);
  }

  friend bool operator!= (eT lhs, hyperdual rhs) {
    return (lhs != rhs.f0);
  }

  friend bool operator!= (hyperdual lhs, eT rhs) {
    return (lhs.f0 != rhs);
  }
};

namespace Eigen {
namespace internal {
    template <typename _Tp>
    struct cast_impl<hyperdual<_Tp>, _Tp>
    {
      static inline _Tp run(const hyperdual<_Tp> & x) {
        return (x.real());
      }
    };
    template <typename _Tp>
    struct cast_impl<hyperdual<_Tp>, int>
    {
      static inline int run(const hyperdual<_Tp> & x) {
        return int(x.real());
      }
    };

}
}

// template <class eT> using hyperducx = hyperdual<std::complex<eT> >;

