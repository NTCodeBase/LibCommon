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

#include <LibCommon/CommonSetup.h>

#define EXPECT_NEAR(x, y, z)

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace Optimization {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
class Problem {
protected:

    bool hasLowerBound_ = false;
    bool hasUpperBound_ = false;

    StdVT<Real_t> lowerBound_;
    StdVT<Real_t> upperBound_;

public:

    Problem() {}

    void setBoxConstraint(StdVT<Real_t> lb, StdVT<Real_t> ub) {
        setLowerBound(lb);
        setUpperBound(ub);
    }

    void setLowerBound(StdVT<Real_t> lb) {
        lowerBound_    = lb;
        hasLowerBound_ = true;
    }

    void setUpperBound(StdVT<Real_t> ub) {
        upperBound_    = ub;
        hasUpperBound_ = true;
    }

    bool hasLowerBound() {
        return hasLowerBound_;
    }

    bool hasUpperBound() {
        return hasUpperBound_;
    }

    StdVT<Real_t> lowerBound() {
        return lowerBound_;
    }

    StdVT<Real_t> upperBound() {
        return upperBound_;
    }

    /**
     * @brief returns objective value in x
     * @details [long description]
     *
     * @param x [description]
     * @return [description]
     */
    virtual Real_t value(const StdVT<Real_t>& x) = 0;

    /**
     * @brief overload value for nice syntax
     * @details [long description]
     *
     * @param x [description]
     * @return [description]
     */
    Real_t operator()(const StdVT<Real_t>& x) {
        return value(x);
    }

    /**
     * @brief returns gradient in x as reference parameter
     * @details should be overwritten by symbolic gradient
     *
     * @param grad [description]
     */
    virtual void gradient(const StdVT<Real_t>& x, StdVT<Real_t>& grad) {
        finiteGradient(x, grad);
    }

    // Compute the value AND the gradient
    virtual Real_t valueGradient(const StdVT<Real_t>& x, StdVT<Real_t>& grad) {
        gradient(x, grad);
        return value(x);
    }

    /**
     * @brief This computes the hessian
     * @details should be overwritten by symbolic hessian, if solver relies on hessian
     */
    //virtual void hessian(const StdVT<Real_t>& x, MatrixX<Real_t>& hessian)
    //{
    //    finiteHessian(x, hessian);
    //}

    virtual bool checkGradient(const StdVT<Real_t>& x, int accuracy = 3) {
        // TODO: check if derived class exists:
        // int(typeid(&Rosenbrock<float>::gradient) == typeid(&Problem<float>::gradient)) == 1 --> overwritten
        const size_t  D = x.size();
        StdVT<Real_t> actual_grad(D);
        StdVT<Real_t> expected_grad(D);
        gradient(x, actual_grad);
        finiteGradient(x, expected_grad, accuracy);

        bool correct = true;

        for(size_t d = 0; d < D; ++d) {
            Real_t scale = std::max((std::max(fabs(actual_grad[d]), fabs(expected_grad[d]))), Real_t(1.));
            EXPECT_NEAR(actual_grad[d], expected_grad[d], 1e-2 * scale);
            if(fabs(actual_grad[d] - expected_grad[d]) > 1e-2 * scale) {
                correct = false;
            }
        }
        return correct;
    }

    //virtual bool checkHessian(const StdVT<Real_t>& x, int accuracy = 3)
    //{
    //    // TODO: check if derived class exists:
    //    // int(typeid(&Rosenbrock<float>::gradient) == typeid(&Problem<float>::gradient)) == 1 --> overwritten
    //    const int D       = x.rows();
    //    bool      correct = true;

    //    MatrixX<Real_t> actual_hessian   = MatrixX<Real_t>::Zero(D, D);
    //    MatrixX<Real_t> expected_hessian = MatrixX<Real_t>::Zero(D, D);
    //    hessian(x, actual_hessian);
    //    finiteHessian(x, expected_hessian, accuracy);
    //    for(int d = 0; d < D; ++d) {
    //        for(int e = 0; e < D; ++e) {
    //            Real_t scale = std::max(static_cast<Real_t>(std::max(fabs(actual_hessian(d, e)), fabs(expected_hessian(d, e)))), (Real_t)1.);
    //            EXPECT_NEAR(actual_hessian(d, e), expected_hessian(d, e), 1e-1 * scale);
    //            if(fabs(actual_hessian(d, e) - expected_hessian(d, e)) > 1e-1 * scale) {
    //                correct = false;
    //            }
    //        }
    //    }
    //    return correct;
    //}

    virtual void finiteGradient(const StdVT<Real_t>& x, StdVT<Real_t>& grad, int accuracy = 0) final {
        // accuracy can be 0, 1, 2, 3
        const Real_t                           eps   = Real_t(2.2204e-6);
        const size_t                           D     = x.size();
        const std::vector<std::vector<Real_t>> coeff =
        { { 1, -1 }, { 1, -8, 8, -1 }, { -1, 9, -45, 45, -9, 1 }, { 3, -32, 168, -672, 672, -168, 32, -3 } };
        const std::vector<std::vector<Real_t>> coeff2 =
        { { 1, -1 }, { -2, -1, 1, 2 }, { -3, -2, -1, 1, 2, 3 }, { -4, -3, -2, -1, 1, 2, 3, 4 } };
        const std::vector<Real_t> dd = { 2, 12, 60, 840 };

        StdVT<Real_t> finiteDiff(D);
        for(size_t d = 0; d < D; d++) {
            finiteDiff[d] = 0;
            for(int s = 0; s < 2 * (accuracy + 1); ++s) {
                StdVT<Real_t> xx = x;
                xx[d]         += coeff2[accuracy][s] * eps;
                finiteDiff[d] += coeff[accuracy][s] * value(xx);
            }
            finiteDiff[d] /= (dd[accuracy] * eps);
        }
        grad = finiteDiff;
    }

    /*
            virtual void finiteHessian(const StdVT<Real_t>& x, MatrixX<Real_t>& hessian, int accuracy = 0) final
            {
                    const Real_t eps = std::numeric_limits<Real_t>::epsilon() * 10e7;
                    const size_t   DIM = x.size();

                    if(accuracy == 0) {
                            for(size_t i = 0; i < DIM; i++) {
                                    for(size_t j = 0; j < DIM; j++) {
                                            StdVT<Real_t> xx = x;
                                            Real_t         f4 = value(xx);
                                            xx[i] += eps;
                                            xx[j] += eps;
                                            Real_t f1 = value(xx);
                                            xx[j] -= eps;
                                            Real_t f2 = value(xx);
                                            xx[j] += eps;
                                            xx[i] -= eps;
                                            Real_t f3 = value(xx);
                                            hessian(i, j) = (f1 - f2 - f3 + f4) / (eps * eps);
                                    }
                            }
                    } else {
                            StdVT<Real_t> xx;
                            for(size_t i = 0; i < DIM; i++) {
                                    for(size_t j = 0; j < DIM; j++) {
                                            Real_t term_1 = 0;
                                            xx = x.eval(); xx[i] += 1 * eps;  xx[j] += -2 * eps;  term_1 += value(xx);
                                            xx = x.eval(); xx[i] += 2 * eps;  xx[j] += -1 * eps;  term_1 += value(xx);
                                            xx = x.eval(); xx[i] += -2 * eps; xx[j] += 1 * eps;   term_1 += value(xx);
                                            xx = x.eval(); xx[i] += -1 * eps; xx[j] += 2 * eps;   term_1 += value(xx);

                                            Real_t term_2 = 0;
                                            xx = x.eval(); xx[i] += -1 * eps; xx[j] += -2 * eps;  term_2 += value(xx);
                                            xx = x.eval(); xx[i] += -2 * eps; xx[j] += -1 * eps;  term_2 += value(xx);
                                            xx = x.eval(); xx[i] += 1 * eps;  xx[j] += 2 * eps;   term_2 += value(xx);
                                            xx = x.eval(); xx[i] += 2 * eps;  xx[j] += 1 * eps;   term_2 += value(xx);

                                            Real_t term_3 = 0;
                                            xx = x.eval(); xx[i] += 2 * eps;  xx[j] += -2 * eps;  term_3 += value(xx);
                                            xx = x.eval(); xx[i] += -2 * eps; xx[j] += 2 * eps;   term_3 += value(xx);
                                            xx = x.eval(); xx[i] += -2 * eps; xx[j] += -2 * eps;  term_3 -= value(xx);
                                            xx = x.eval(); xx[i] += 2 * eps;  xx[j] += 2 * eps;   term_3 -= value(xx);

                                            Real_t term_4 = 0;
                                            xx = x.eval(); xx[i] += -1 * eps; xx[j] += -1 * eps;  term_4 += value(xx);
                                            xx = x.eval(); xx[i] += 1 * eps;  xx[j] += 1 * eps;   term_4 += value(xx);
                                            xx = x.eval(); xx[i] += 1 * eps;  xx[j] += -1 * eps;  term_4 -= value(xx);
                                            xx = x.eval(); xx[i] += -1 * eps; xx[j] += 1 * eps;   term_4 -= value(xx);

                                            hessian(i, j) = (-63 * term_1 + 63 * term_2 + 44 * term_3 + 74 * term_4) / (600.0 * eps * eps);
                                    }
                            }
                    }
            }*/
};
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace Optimization
