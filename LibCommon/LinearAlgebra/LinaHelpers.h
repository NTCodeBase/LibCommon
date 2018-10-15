//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
//    .--------------------------------------------------.
//    |  This file is part of NTGraphics                 |
//    |  Created 2018 by NT (https://ttnghia.github.io)  |
//    '--------------------------------------------------'
//                            \o/
//                             |
//                            / |
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#pragma once

#include <random>

#include <LibCommon/CommonSetup.h>
#include <LibCommon/Math/MathHelpers.h>
#include <LibCommon/Utils/NumberHelpers.h>
#include <LibCommon/LinearAlgebra/ImplicitQRSVD.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace LinaHelpers {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
bool hasValidElements(const VecX<N, Real_t>& vec) {
    for(Int i = 0; i < N; ++i) {
        if(!NumberHelpers::isValidNumber(vec[i])) {
            return false;
        }
    }
    return true;
}

template<Int N, class Real_t>
bool hasValidElements(const MatXxX<N, Real_t>& mat) {
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            if(!NumberHelpers::isValidNumber(mat[i][j])) {
                return false;
            }
        }
    }
    return true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
Real_t maxAbs(const MatXxX<N, Real_t>& mat) {
    Real_t result = Real_t(0);
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            result = MathHelpers::max(result, fabs(mat[i][j]));
        }
    }
    return result;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
Real_t norm2(const MatXxX<N, Real_t>& mat) {
    Real_t prod = Real_t(0);
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            prod += mat[i][j] * mat[i][j];
        }
    }
    return sqrt(prod);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T, class S>
void fill(MatXxX<N, T>& mat, S x) {
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            mat[i][j] = T(x);
        }
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
Real_t trace(const MatXxX<N, Real_t>& mat) {
    Real_t prod = Real_t(0);
    for(Int i = 0; i < N; ++i) {
        prod += mat[i][i];
    }
    return prod;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
MatXxX<N, Real_t> dev(const MatXxX<N, Real_t>& mat) {
    return mat - MatXxX<N, Real_t>(LinaHelpers::trace<Real_t>(mat) / Real_t(N));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void sumToDiag(MatXxX<N, Real_t>& mat, Real_t c) {
    for(Int i = 0; i < N; ++i) {
        mat[i][i] += c;
    }
}

template<Int N, class Real_t>
void sum1ToDiag(MatXxX<N, Real_t>& mat) {
    for(Int i = 0; i < N; ++i) {
        mat[i][i] += Real_t(1.0);
    }
}

template<Int N, class Real_t>
MatXxX<N, Real_t> getDiagSum(const MatXxX<N, Real_t>& mat, Real_t c) {
    auto result = mat;
    for(Int i = 0; i < N; ++i) {
        result[i][i] += c;
    }
    return result;
}

template<Int N, class Real_t>
MatXxX<N, Real_t> getDiagSum(const MatXxX<N, Real_t>& mat, const VecX<N, Real_t>& c) {
    auto result = mat;
    for(Int i = 0; i < N; ++i) {
        result[i][i] += c[i];
    }
    return result;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
VecX<N, Real_t> extractDiag(const MatXxX<N, Real_t>& mat) {
    VecX<N, Real_t> diag;
    for(Int i = 0; i < N; ++i) {
        diag[i] = mat[i][i];
    }
    return diag;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
MatXxX<N, Real_t> diagMatrix(const VecX<N, Real_t>& diag) {
    MatXxX<N, Real_t> mat(0);
    for(Int i = 0; i < N; ++i) {
        mat[i][i] = diag[i];
    }
    return mat;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void diagProduct(MatXxX<N, Real_t>& mat, const VecX<N, Real_t>& vec) {
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            mat[i][j] *= vec[i];
        }
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//Matrix * Matrix^-1
template<Int N, class Real_t>
void diagProductInv(MatXxX<N, Real_t>& mat, const VecX<N, Real_t>& vec) {
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            mat[i][j] /= vec[i];
        }
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
Real_t frobeniusInnerProduct(const MatXxX<N, Real_t>& m1, const MatXxX<N, Real_t>& m2) {
    Real_t prod = Real_t(0);
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            prod += m1[i][j] * m2[i][j];
        }
    }
    return prod;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
VecX<N, Real_t> innerProduct(const VecX<N, Real_t>& vec, const MatXxX<N, Real_t>& mat) {
    VecX<N, Real_t> prod(0);
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            prod[i] += vec[j] * mat[j][i];
        }
    }
    return prod;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
MatXxX<N, Real_t> innerProduct(const MatXxX<N, Real_t>& m1, const MatXxX<N, Real_t>& m2) {
    MatXxX<N, Real_t> prod(0);
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            for(Int k = 0; k < N; ++k) {
                prod[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }
    return prod;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// TODO: check row-col major
template<class Real_t>
Mat2x2<Real_t> cofactor(const Mat2x2<Real_t>& mat) {
    return Mat2x2<Real_t>(mat[1][1], -mat[0][1],
                          -mat[1][0], mat[0][0]);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// TODO: test value row-col major
template<class Real_t>
Real_t elementCofactor(const Mat3x3<Real_t>& mat, Int x, Int y) {
    Real_t         cofactor_v;
    Real_t         minor;
    Mat2x2<Real_t> minor_mat;

    minor_mat[0][0] = mat[(x + 1) % 3][(y + 1) % 3];
    minor_mat[1][1] = mat[(x + 2) % 3][(y + 2) % 3];
    minor_mat[0][1] = mat[(x + 1) % 3][(y + 2) % 3];
    minor_mat[1][0] = mat[(x + 2) % 3][(y + 1) % 3];
    minor           = glm::determinant(minor_mat);

    cofactor_v = ((x + y) & 1) ? -minor : minor;
    if(y == 1) {
        cofactor_v = -cofactor_v;
    }
    return cofactor_v;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
Mat3x3<Real_t> cofactor(const Mat3x3<Real_t>& mat) {
    Mat2x2<Real_t> result;
    for(Int i = 0; i < 2; ++i) {
        for(Int j = 0; j < 2; ++j) {
            result[i][j] = elementCofactor(mat, i, j);
        }
    }
    return result;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
Real_t vonMisesPlaneStress(const Mat3x3<Real_t>& mat) {
    const Real_t vm = mat[0][0] * mat[0][0] + mat[1][1] * mat[1][1] + mat[2][2] * mat[2][2] -
                      mat[0][0] * mat[1][1] - mat[1][1] * mat[2][2] - mat[2][2] * mat[0][0] +
                      (mat[0][1] * mat[1][0] + mat[1][2] * mat[2][1] + mat[2][0] * mat[0][2]) * Real_t(3.0);

    return vm > 0 ? std::sqrt(vm) : 0;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
MatXxX<N, Real_t> randMatrix(Real_t minVal = Real_t(0), Real_t maxVal = Real_t(1.0)) {
    std::random_device                     rd;
    std::mt19937                           gen(rd());
    std::uniform_real_distribution<Real_t> dis(minVal, maxVal);

    MatXxX<N, Real_t> result;
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            result[i][j] = dis(gen);
        }
    }
    return result;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class SizeType, class Real_t>
StdVT<MatXxX<N, Real_t>> randVecMatrices(SizeType size, Real_t minVal = Real_t(0), Real_t maxVal = Real_t(1.0)) {
    std::random_device                     rd;
    std::mt19937                           gen(rd());
    std::uniform_real_distribution<Real_t> dis(minVal, maxVal);

    StdVT<MatXxX<N, Real_t>> results;
    results.resize(size);
    for(SizeType idx = 0; idx < size; ++idx) {
        for(Int i = 0; i < N; ++i) {
            for(Int j = 0; j < N; ++j) {
                results[idx][i][j] = dis(gen);
            }
        }
    }
    return results;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
auto orientedSVD(const MatXxX<N, Real_t>& M) {
    MatXxX<N, Real_t> U, Vt;
    VecX<N, Real_t>   S;

    QRSVD::svd(M, U, S, Vt);
    Vt = glm::transpose(Vt);

    MatXxX<N, Real_t> J = MatXxX<N, Real_t>(1.0);
    J[N - 1][N - 1] = Real_t(-1.0);

    // Check for inversion
    if(glm::determinant(U) < Real_t(0)) {
        U         = U * J;
        S[N - 1] *= Real_t(-1.0);
    }
    if(glm::determinant(Vt) < Real_t(0)) {
        Vt        = J * Vt;
        S[N - 1] *= Real_t(-1.0);
    }

    return std::make_tuple(U, S, Vt);
} // end oriented svd

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
void QRDifferential(const Mat2x2<Real_t>& Q, const Mat2x2<Real_t>& R, const Mat2x2<Real_t>& dF, Mat2x2<Real_t>& dQ, Mat2x2<Real_t>& dR) {
    __NT_REQUIRE(R[0][0] != 0);
    Mat2x2<Real_t> QtdF = glm::transpose(Q) * dF;
    Real_t         a    = -QtdF[0][1] / R[0][0];
    Mat2x2<Real_t> QtdQ(0, -a, a, 0);
    dQ = Q * QtdQ;
    dR = glm::transpose(Q) * dF - QtdQ * R;
}

template<class Real_t>
void QRDifferential(const Mat3x3<Real_t>& Q, const Mat3x3<Real_t>& R, const Mat3x3<Real_t>& dF, Mat3x3<Real_t>& dQ, Mat3x3<Real_t>& dR) {
    __NT_REQUIRE(R[0][0] != 0 && R[1][1] != 0);
    Mat3x3<Real_t> QtdF = glm::transpose(Q) * dF;
    Real_t         w3   = QtdF[0][1] / R[0][0];
    Real_t         w2   = -QtdF[0][2] / R[0][0];
    Real_t         w1   = (QtdF[1][2] + w2 * R[1][0]) / R[1][1];
    Mat3x3<Real_t> QtdQ(0, w3, -w2, -w3, 0, w1, w2, -w1, 0);
    dQ = Q * QtdQ;
    dR = glm::transpose(Q) * dF - QtdQ * R;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
auto symmetryDecomposition(const MatXxX<N, Real_t>& M) {
    MatXxX<N, Real_t> symComp, skewSymComp;
    auto              Mt = glm::transpose(M);

    symComp     = Real_t(0.5) * (M + Mt);
    skewSymComp = Real_t(0.5) * (M - Mt);
    return std::make_tuple(symComp, skewSymComp);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
auto extractFiberCotangentStress(const MatXxX<N, Real_t>& VP, const MatXxX<N, Real_t>& F) {
    if constexpr (N == 2) {
        return glm::outerProduct(VP[1], F[1]);
    } else {
        return MatMxN<3, 2, Real_t>(VP[1], VP[2]) * glm::transpose(MatMxN<3, 2, Real_t>(F[1], F[2]));
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
MatXxX<N, Real_t> getOrthogonalSystem(const VecX<N, Real_t>& d1) {
    if constexpr (N == 2) {
        MatXxX<N, Real_t> M;
        M[0] = glm::normalize(d1);
        M[1] = glm::normalize(Vec2<Real_t>(-d1.y, d1.x));
        return M;
    } else {
        const Vec3<Real_t> basis[3] = { Vec3<Real_t>(1, 0, 0), Vec3<Real_t>(0, 1, 0), Vec3<Real_t>(0, 0, 1) };
        Vec3<Real_t>       u;
        MatXxX<N, Real_t>  M;
        M[0] = glm::normalize(d1);

        for(UInt i = 0; i < 3; ++i) {
            if(auto de = glm::dot(M[0], basis[i]); de > Real_t(0.1)) {
                u = basis[i];
                break;
            }
        }
        auto v = glm::cross(d1, u);
        u = glm::cross(v, d1);

        M[1] = glm::normalize(u);
        M[2] = glm::normalize(v);
        return M;
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Compute the rotation quaternion to rotate from position p1 to p2
template<class T>
auto angleAxisFromPositions(const Vec3<T>& p1, const Vec3<T>& p2) {
    auto tmp = glm::dot(glm::normalize(p1), glm::normalize(p2));
    if(tmp > T(1)) { tmp = T(1); }
    auto angle = std::acos(tmp);
    auto axis  = glm::normalize(glm::cross(p1, p2));
    return std::make_pair(angle, axis);
}

template<class T>
Quat<T> quaternionFromPositions(const Vec3<T>& p1, const Vec3<T>& p2) {
    auto [angle, axis] = angleAxisFromPositions(p1, p2);
    return glm::angleAxis(angle, axis);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Quaternion>
auto quaternionMatrix(const Quaternion& q) {
    return Mat4x4<Quaternion::value_type>(q.w, q.z, -q.y, -q.x, // col 0
                                          -q.z, q.w, q.x, -q.y, // col 1
                                          q.y, -q.x, q.w, -q.z, // col 2
                                          q.x, q.y, q.z, q.w);  // col 3
}

template<class Quaternion>
auto quaternionMatrixTransposed3x4(const Quaternion& q) {
    return MatMxN<3, 4, Quaternion::value_type>(q.w, -q.z, q.y,    // col 0
                                                q.z, q.w, -q.x,    // col 1
                                                -q.y, q.x, q.w,    // col 2
                                                -q.x, -q.y, -q.z); // col 3
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T> Quat<T> vecXToQuat(const Vec2<T>& q) { return Quat<T>(q.y, 0, 0, q.x); }
template<class T> Quat<T> vecXToQuat(const Vec4<T>& q) { return Quat<T>(q.w, q.x, q.y, q.z); }
template<Int N, class T> VecX<N, T> quatToVecX(const Quat<T>& q) { if constexpr (N == 2) { return Vec2<T>(q.z, q.w); } else { return Vec4<T>(q.x, q.y, q.z, q.w); } }

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace LinaHelpers
