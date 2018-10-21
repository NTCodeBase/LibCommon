//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
//    .--------------------------------------------------.
//    |  This file is part of NTCodeBase                 |
//    |  Created 2018 by NT (https://ttnghia.github.io)  |
//    '--------------------------------------------------'
//                            \o/
//                             |
//                            / |
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#pragma once

#include <LibCommon/ParallelHelpers/ParallelBLAS.h>
#include <LibCommon/ParallelHelpers/ParallelSTL.h>
#include <LibCommon/Math/Optimization/MoreThuente.h>
#include <LibCommon/Math/Optimization/Problem.h>

#include <iostream>
#include <algorithm>
#include <cmath>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace Banana::Optimization {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
class MatrixX {
public:
    template<class IndexType>
    MatrixX(IndexType nRows, IndexType nCols, Real_t defaultVal = Real_t(0)) : m_nRows(static_cast<UInt>(nRows)), m_nCols(static_cast<UInt>(nCols)) {
        m_Data.resize(nCols);
        for(auto& col : m_Data) {
            col.resize(m_nRows, defaultVal);
        }
    }

    auto nRows() const { return m_nRows; }
    auto nCols() const { return m_nCols; }

    template<class IndexType> const auto& col(IndexType col) const { return m_Data[col]; }
    template<class IndexType> auto& col(IndexType col) { return m_Data[col]; }

    const auto& data() const { return m_Data; }
    auto& data() { return m_Data; }

private:
    UInt                 m_nRows = 0;
    UInt                 m_nCols = 0;
    StdVT<StdVT<Real_t>> m_Data;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<typename Real_t, Int Order>
class ISolver {
public:
    ISolver() = default;

    /**
     * @brief minimize an objective function given a gradient
     * @details this is just the abstract interface
     *
     * @param x0 starting point
     * @param funObjective objective function
     */
    virtual void minimize(Problem<Real_t>& objFunc, StdVT<Real_t>& x0) = 0;
    ////////////////////////////////////////////////////////////////////////////////
    auto& gradTolerance() { return m_GradTol; }
    auto& initHessian() { return m_InitHess; }
    auto& maxIter() { return m_MaxIter; }
    ////////////////////////////////////////////////////////////////////////////////
    auto        nIters() const { return m_nIters; }
    const auto& objVals() const { return m_ObjVals; }
    const auto& runtimes()  const { return m_Runtimes; }
protected:
    size_t        m_nIters;
    StdVT<Real_t> m_ObjVals;
    StdVT<Real_t> m_Runtimes;

    Real_t m_GradTol  = Real_t(1e-4);
    Real_t m_InitHess = Real_t(1.0); // only used by lbfgs
    size_t m_MaxIter  = 100000;
};
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

/**
 * @brief  LBFGS implementation based on Nocedal & Wright Numerical Optimization book (Section 7.2)
 * @tparam T scalar type
 * @tparam P problem type
 * @tparam Order order of solver
 */

template<class Real_t>
class LBFGSSolver : public ISolver<Real_t, 1> {
public:
    void minimize(Problem<Real_t>& objFunc, StdVT<Real_t>& x0) {
        size_t m     = std::min(this->m_MaxIter, size_t(10));
        size_t nVars = static_cast<UInt>(x0.size());
        Real_t eps_g = this->m_GradTol;
        Real_t eps_x = Real_t(1e-8);

        auto s = MatrixX<Real_t>(nVars, m);
        auto y = MatrixX<Real_t>(nVars, m);

        StdVT<Real_t> alpha = StdVT<Real_t>(m, 0);
        StdVT<Real_t> rho   = StdVT<Real_t>(m, 0);
        StdVT<Real_t> grad(nVars), q(nVars), grad_old(nVars), x_old(nVars);

        //	Real_t f = objFunc.value(x0);
        Real_t f              = objFunc.valueGradient(x0, grad);
        Real_t gamma_k        = this->m_InitHess;
        Real_t gradNorm       = 0;
        Real_t alpha_init     = std::min(Real_t(1.0), Real_t(1.0) / ParallelSTL::maxAbs(grad));
        size_t globIter       = 0;
        size_t maxiter        = this->m_MaxIter;
        Real_t new_hess_guess = 1.0; // only changed if we converged to a solution

        for(size_t k = 0; k < maxiter; k++) {
            x_old    = x0;
            grad_old = grad;
            q        = grad;
            globIter++;

            //L - BFGS first - loop recursion
            size_t iter = std::min(m, k);
            for(Int i = Int(iter) - 1; i >= 0; --i) {
                rho[i]   = Real_t(1.0) / ParallelBLAS::dotProduct(s.col(i), y.col(i));
                alpha[i] = rho[i] * ParallelBLAS::dotProduct(s.col(i), q);
                ParallelBLAS::addScaled(-alpha[i], y.col(i), q);
                //q = q - alpha[i] * y.col(i);
            }

            //L - BFGS second - loop recursion
            //q = gamma_k * q;
            ParallelBLAS::scale(gamma_k, q);
            for(size_t i = 0; i < iter; ++i) {
                Real_t beta = rho[i] * ParallelBLAS::dotProduct(q, y.col(i));
                ParallelBLAS::addScaled((alpha[i] - beta), s.col(i), q);
                //q = q + (alpha[i] - beta) * s.col(i);
            }

            // is there a descent
            Real_t dir = ParallelBLAS::dotProduct(q, grad);
            if(dir < Real_t(1e-4)) {
                q          = grad;
                maxiter   -= k;
                k          = 0;
                alpha_init = std::min(Real_t(1.0), Real_t(1.0) / ParallelSTL::maxAbs(grad));
            }

            const Real_t rate = MoreThuente<Real_t, decltype(objFunc), 1>::linesearch(x0, ParallelBLAS::multiply(Real_t(-1.0), q), objFunc, alpha_init);
            //		const Real_t rate = linesearch(objFunc, x0, -q, f, grad, 1.0);

            ParallelBLAS::addScaled(-rate, q, x0);
            //x0 = x0 - rate * q;

            if(ParallelBLAS::norm2(ParallelBLAS::minus(x_old, x0)) < eps_x) {
                //			std::cout << "x diff norm: " << (x_old - x0).squaredNorm() << std::endl;
                break;
            } // usually this is a problem so exit

            //		f = objFunc.value(x0);
            f = objFunc.valueGradient(x0, grad);

            gradNorm = ParallelSTL::maxAbs(grad);
            if(gradNorm < eps_g) {
                // Only change hessian guess if we break out the loop via convergence.
                //			std::cout << "grad norm: " << gradNorm << std::endl;
                new_hess_guess = gamma_k;
                break;
            }

            StdVT<Real_t> s_temp = ParallelBLAS::minus(x0, x_old);
            StdVT<Real_t> y_temp = ParallelBLAS::minus(grad, grad_old);
            //StdVT<Real_t> s_temp = x0 - x_old;
            //StdVT<Real_t> y_temp = grad - grad_old;

            // update the history
            if(k < m) {
                s.col(k) = s_temp;
                y.col(k) = y_temp;
            } else {
                std::rotate(s.data().begin(), s.data().begin() + 1, s.data().end());
                s.col(s.nCols() - 1) = s_temp;

                std::rotate(y.data().begin(), y.data().begin() + 1, y.data().end());
                y.col(y.nCols() - 1) = y_temp;

                /*s.leftCols(_m - 1) = s.rightCols(_m - 1).eval();
                   s.rightCols(1)     = s_temp;
                   y.leftCols(_m - 1) = y.rightCols(_m - 1).eval();
                   y.rightCols(1)     = y_temp;*/
            }

            gamma_k    = ParallelBLAS::dotProduct(s_temp, y_temp) / ParallelBLAS::dotProduct(y_temp, y_temp);
            alpha_init = 1.0;
        }

        this->m_nIters   = globIter;
        this->m_InitHess = new_hess_guess;
    }
};
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace Banana::Optimization
