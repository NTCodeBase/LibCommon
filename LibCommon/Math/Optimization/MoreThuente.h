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
#include <cmath>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace Optimization {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t, typename P, Int Order>
class MoreThuente {
public:

    /**
     * @brief use MoreThuente Rule for (strong) Wolfe conditiions
     * @details [long description]
     *
     * @param searchDir search direction for next update step
     * @param objFunc handle to problem
     *
     * @return step-width
     */

    static Real_t linesearch(const StdVT<Real_t>& x, const StdVT<Real_t>& searchDir, P& objFunc, const Real_t alpha_init = Real_t(1.0)) {
        // assume step width
        Real_t        ak   = alpha_init;
        StdVT<Real_t> g    = x;
        Real_t        fval = objFunc.valueGradient(x, g);
        StdVT<Real_t> xx   = x;
        StdVT<Real_t> s    = searchDir;
        ////////////////////////////////////////////////////////////////////////////////
        cvsrch(objFunc, xx, fval, g, ak, s);
        return ak;
    }

    static Int cvsrch(P& objFunc, StdVT<Real_t>& x, Real_t f, StdVT<Real_t>& g, Real_t& stp, StdVT<Real_t>& s) {
        Int          info   = 0;
        Int          infoc  = 1;
        const Real_t xtol   = Real_t(1e-15);
        const Real_t ftol   = Real_t(1e-4);
        const Real_t gtol   = Real_t(1e-2);
        const Real_t stpmin = Real_t(1e-15);
        const Real_t stpmax = Real_t(1e15);
        const Real_t xtrapf = Real_t(4);
        const Int    maxfev = 20;
        Int          nfev   = 0;

        Real_t dginit = ParallelBLAS::dotProduct(g, s);
        if(dginit >= 0) {
            // no descent direction
            // TODO: handle this case
            return -1;
        }

        bool brackt = false;
        bool stage1 = true;

        Real_t        finit  = f;
        Real_t        dgtest = ftol * dginit;
        Real_t        width  = stpmax - stpmin;
        Real_t        width1 = 2 * width;
        StdVT<Real_t> wa     = x;

        Real_t stx = Real_t(0.0);
        Real_t fx  = finit;
        Real_t dgx = dginit;
        Real_t sty = Real_t(0.0);
        Real_t fy  = finit;
        Real_t dgy = dginit;

        Real_t stmin;
        Real_t stmax;

        while(true) {
            // make sure we stay in the interval when setting min/max-step-width
            if(brackt) {
                stmin = std::min(stx, sty);
                stmax = std::max(stx, sty);
            } else {
                stmin = stx;
                stmax = stp + xtrapf * (stp - stx);
            }

            // Force the step to be within the bounds stpmax and stpmin.
            stp = std::max(stp, stpmin);
            stp = std::min(stp, stpmax);

            // Oops, let us return the last reliable values
            if((brackt && (stp <= stmin || stp >= stmax)) |
               (nfev >= maxfev - 1) |
               (infoc == 0) | (brackt & (stmax - stmin <= xtol * stmax))) {
                stp = stx;
            }

            // test new point
            x = wa;
            ParallelBLAS::addScaled(stp, s, x);
            //x = wa + stp * s;

            //    f = objFunc.value(x);
            //    objFunc.gradient(x, g);
            f = objFunc.valueGradient(x, g);

            nfev++;
            Real_t dg     = ParallelBLAS::dotProduct(g, s);
            Real_t ftest1 = finit + stp * dgtest;

            // all possible convergence tests
            if((brackt & ((stp <= stmin) | (stp >= stmax))) | (infoc == 0)) {
                info = 6;
            }

            if((stp == stpmax) & (f <= ftest1) & (dg <= dgtest)) {
                info = 5;
            }

            if((stp == stpmin) & ((f > ftest1) | (dg >= dgtest))) {
                info = 4;
            }

            if(nfev >= maxfev) {
                info = 3;
            }

            if(brackt & (stmax - stmin <= xtol * stmax)) {
                info = 2;
            }

            if((f <= ftest1) & (fabs(dg) <= gtol * (-dginit))) {
                info = 1;
            }

            // terminate when convergence reached
            if(info != 0) {
                return -1;
            }

            if(stage1 & (f <= ftest1) & (dg >= std::min(ftol, gtol) * dginit)) {
                stage1 = false;
            }

            if(stage1 & (f <= fx) & (f > ftest1)) {
                Real_t fm   = f - stp * dgtest;
                Real_t fxm  = fx - stx * dgtest;
                Real_t fym  = fy - sty * dgtest;
                Real_t dgm  = dg - dgtest;
                Real_t dgxm = dgx - dgtest;
                Real_t dgym = dgy - dgtest;

                cstep(stx, fxm, dgxm, sty, fym, dgym, stp, fm, dgm, brackt, stmin, stmax, infoc);

                fx  = fxm + stx * dgtest;
                fy  = fym + sty * dgtest;
                dgx = dgxm + dgtest;
                dgy = dgym + dgtest;
            } else {
                // this is ugly and some variables should be moved to the class scope
                cstep(stx, fx, dgx, sty, fy, dgy, stp, f, dg, brackt, stmin, stmax, infoc);
            }

            if(brackt) {
                if(fabs(sty - stx) >= Real_t(0.66) * width1) {
                    stp = stx + Real_t(0.5) * (sty - stx);
                }
                width1 = width;
                width  = fabs(sty - stx);
            }
        }

        return 0;
    }

    static Int cstep(Real_t& stx, Real_t& fx, Real_t& dx, Real_t& sty, Real_t& fy, Real_t& dy, Real_t& stp,
                     Real_t& fp, Real_t& dp, bool& brackt, Real_t& stpmin, Real_t& stpmax, Int& info) {
        info = 0;
        bool bound = false;

        // Check the input parameters for errors.
        if((brackt & ((stp <= std::min(stx, sty)) | (stp >= std::max(stx, sty)))) | (dx * (stp - stx) >= 0.0)
           | (stpmax < stpmin)) {
            return -1;
        }

        Real_t sgnd = dp * (dx / fabs(dx));

        Real_t stpf = 0;
        Real_t stpc = 0;
        Real_t stpq = 0;

        if(fp > fx) {
            info  = 1;
            bound = true;
            Real_t theta = Real_t(3.0) * (fx - fp) / (stp - stx) + dx + dp;
            Real_t s     = std::max(theta, std::max(dx, dp));
            Real_t gamma = s * sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
            if(stp < stx) {
                gamma = -gamma;
            }
            Real_t p = (gamma - dx) + theta;
            Real_t q = ((gamma - dx) + gamma) + dp;
            Real_t r = p / q;
            stpc = stx + r * (stp - stx);
            stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / Real_t(2.0)) * (stp - stx);
            if(fabs(stpc - stx) < fabs(stpq - stx)) {
                stpf = stpc;
            } else {
                stpf = stpc + (stpq - stpc) / 2;
            }
            brackt = true;
        } else if(sgnd < 0.0) {
            info  = 2;
            bound = false;
            Real_t theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
            Real_t s     = std::max(theta, std::max(dx, dp));
            Real_t gamma = s * sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
            if(stp > stx) {
                gamma = -gamma;
            }

            Real_t p = (gamma - dp) + theta;
            Real_t q = ((gamma - dp) + gamma) + dx;
            Real_t r = p / q;
            stpc = stp + r * (stx - stp);
            stpq = stp + (dp / (dp - dx)) * (stx - stp);
            if(fabs(stpc - stp) > fabs(stpq - stp)) {
                stpf = stpc;
            } else {
                stpf = stpq;
            }
            brackt = true;
        } else if(fabs(dp) < fabs(dx)) {
            info  = 3;
            bound = 1;
            Real_t theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
            Real_t s     = std::max(theta, std::max(dx, dp));
            Real_t gamma = s * sqrt(std::max(static_cast<Real_t>(0.), (theta / s) * (theta / s) - (dx / s) * (dp / s)));
            if(stp > stx) {
                gamma = -gamma;
            }
            Real_t p = (gamma - dp) + theta;
            Real_t q = (gamma + (dx - dp)) + gamma;
            Real_t r = p / q;
            if((r < 0.0) & (gamma != 0.0)) {
                stpc = stp + r * (stx - stp);
            } else if(stp > stx) {
                stpc = stpmax;
            } else {
                stpc = stpmin;
            }
            stpq = stp + (dp / (dp - dx)) * (stx - stp);
            if(brackt) {
                if(fabs(stp - stpc) < fabs(stp - stpq)) {
                    stpf = stpc;
                } else {
                    stpf = stpq;
                }
            } else {
                if(fabs(stp - stpc) > fabs(stp - stpq)) {
                    stpf = stpc;
                } else {
                    stpf = stpq;
                }
            }
        } else {
            info  = 4;
            bound = false;
            if(brackt) {
                Real_t theta = 3 * (fp - fy) / (sty - stp) + dy + dp;
                Real_t s     = std::max(theta, std::max(dy, dp));
                Real_t gamma = s * sqrt((theta / s) * (theta / s) - (dy / s) * (dp / s));
                if(stp > sty) {
                    gamma = -gamma;
                }

                Real_t p = (gamma - dp) + theta;
                Real_t q = ((gamma - dp) + gamma) + dy;
                Real_t r = p / q;
                stpc = stp + r * (sty - stp);
                stpf = stpc;
            } else if(stp > stx) {
                stpf = stpmax;
            } else {
                stpf = stpmin;
            }
        }

        if(fp > fx) {
            sty = stp;
            fy  = fp;
            dy  = dp;
        } else {
            if(sgnd < 0.0) {
                sty = stx;
                fy  = fx;
                dy  = dx;
            }

            stx = stp;
            fx  = fp;
            dx  = dp;
        }

        stpf = std::min(stpmax, stpf);
        stpf = std::max(stpmin, stpf);
        stp  = stpf;

        if(brackt & bound) {
            if(sty > stx) {
                stp = std::min(stx + static_cast<Real_t>(0.66) * (sty - stx), stp);
            } else {
                stp = std::max(stx + static_cast<Real_t>(0.66) * (sty - stx), stp);
            }
        }

        return 0;
    }
};
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace Optimization
