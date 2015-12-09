#ifndef AUTODIFF_H
#define AUTODIFF_H

#include <Eigen/Dense>
#ifndef MATLAB
#include "../gtest/gtest.h"
#else
#define EXPECT_NEAR(x, y, z)
#endif /* MATLAB */

#include "meta.h"

namespace cppoptlib {

template<template<typename> class P, typename TT>
class AutoDiff : public Problem<TT> {
 private:

  P<hyperdual<TT> > *cascadeProblem;
  P<TT >            *originalProblem;

 public:

  AutoDiff() {
    cascadeProblem  = new P<hyperdual<TT> >();
    originalProblem = new P<TT >();
  }

  TT value(const Vector<TT> &x) {
    return originalProblem->value(x);
  }

  void gradient(const Vector<TT> &x, Vector<TT> &grad) {

    const int D = x.rows();
    hyperdual<TT> ans;

    Vector<hyperdual<TT> > xx(D);
    for (int i = 0; i < D; ++i)
    {
      xx[i] = x[i];
    }
    for (int i = 0; i < D; ++i) {
      xx[i].wrt();
      ans = cascadeProblem->value(xx);
      grad[i] = ans.grad();
      xx[i].nowrt();

    }

  }

  void hessian(const Vector<TT> &x, Matrix<TT> &hes) {

    const int D = x.rows();
    hyperdual<TT> ans;

    Vector<hyperdual<TT> > xx(D);
    for (int i = 0; i < D; ++i)
    {
      xx[i] = x[i];

    }

    for (int i = 0; i < D; ++i)
    {
      // xx[i].wrt(3);
      for (int j = 0; j < D; ++j)
      {
        xx[i].wrt(1);
        // xx[i].wrt(3);
        // xx[i].wrt(2);
        // xx[j].wrt(1);
        xx[j].wrt(2);
        xx[j].wrt(3);
        ans = cascadeProblem->value(xx);
        hes(i, j) = ans.hessian();
        for (int k = 0; k < D; ++k)
        {
          xx[k].nowrt();
        }
      }
    }
  }

};

}

#endif /* AUTODIFF_H */
