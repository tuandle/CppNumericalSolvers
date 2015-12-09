#include <iostream>
#include "../../include/cppoptlib/meta.h"
#include "../../include/cppoptlib/problem.h"
#include "../../include/cppoptlib/hyperdual.h"
#include "../../include/cppoptlib/autodiff.h"
#include "../../include/cppoptlib/solver/bfgssolver.h"

// nolintnextline
using namespace cppoptlib;

template<typename T>
class Rosenbrock : public Problem<T> {
  public:
    // this is just the objective (NOT optional)
    T value(const Vector<T> &x) {
        const T t1 = (1 - x[0]);
        const T t2 = (x[1] - x[0] * x[0]);
        return   t1 * t1 + 100 * t2 * t2;
    }

    void gradient(const Vector<T> &x, Vector<T> &grad) {
        grad[0]  = -2 * (1 - x[0]) + 200 * (x[1] - x[0] * x[0]) * (-2 * x[0]);
        grad[1]  =                   200 * (x[1] - x[0] * x[0]);
    }

    void hessian(const Vector<T> &x, Matrix<T> & hessian) {
        hessian(0, 0) = 1200 * x[0] * x[0] - 400 * x[1] + 1;
        hessian(0, 1) = -400 * x[0];
        hessian(1, 0) = -400 * x[0];
        hessian(1, 1) = 200;
    }
};

template<typename T>
class RosenbrockObjective : public Problem<T> {
  public:
    // this is just the objective (NOT optional)
    T value(const Vector<T> &x) {
        const T t1 = (1 - x[0]);
        const T t2 = (x[1] - x[0] * x[0]);
        return   t1 * t1 + 100 * t2 * t2;
    }
};

int main(int argc, char const *argv[]) {

    typedef double T;

    // we diff our function in x0
    Vector<T> x0(2); x0 << -1, 2;

    // symbolic user-given gradient
    Rosenbrock<T> f_sym;
    Vector<T> sym_grad(2);    sym_grad << 0, 0;
    Matrix<T> sym_hess(2, 2); sym_hess << 0, 0, 0, 0;
    f_sym.gradient(x0, sym_grad);
    f_sym.hessian(x0, sym_hess);
    std::cout << sym_grad << std::endl;
    std::cout << sym_hess << std::endl;

    // use auto-diff
    AutoDiff<RosenbrockObjective, T> f_auto;
    Vector<T> auto_grad(2);    auto_grad << 0, 0;
    Matrix<T> auto_hess(2, 2); auto_hess << 0, 0, 0, 0;
    f_auto.gradient(x0, auto_grad);
    f_auto.hessian(x0, auto_hess);
    std::cout << auto_grad << std::endl;
    std::cout << auto_hess << std::endl;

    return 0;
}
