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

#include <random>

#include <LibCommon/CommonSetup.h>
#include <LibCommon/Math/MathHelpers.h>
#include <LibCommon/Utils/NumberHelpers.h>
#include <LibCommon/LinearAlgebra/ImplicitQRSVD.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NTCodeBase::LinaHelpers {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
bool hasValidElements(const VecX<N, T>& vec) {
    for(Int i = 0; i < N; ++i) {
        if(!NumberHelpers::isValidNumber(vec[i])) {
            return false;
        }
    }
    return true;
}

template<Int N, class T>
bool hasValidElements(const MatXxX<N, T>& mat) {
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
// return 0 if ||v|| is too small
template<Int N, class T>
VecX<N, T> safeNormalize(const VecX<N, T>& v) {
    T l2 = glm::length2(v);
    if(l2 > T(1e-20)) {
        return v / std::sqrt(l2);
    }
    return VecX<N, T>(0);
}

// return v0 if ||v|| is too small
template<Int N, class T>
VecX<N, T> safeNormalize(const VecX<N, T>& v, const VecX<N, T>& v0) {
    T l2 = glm::length2(v);
    if(l2 > T(1e-20)) {
        return v / std::sqrt(l2);
    }
    return v0;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
T maxAbs(const MatXxX<N, T>& mat) {
    T result = T(0);
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            result = MathHelpers::max(result, std::abs(mat[i][j]));
        }
    }
    return result;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
T norm2(const MatXxX<N, T>& mat) {
    T prod = T(0);
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            prod += mat[i][j] * mat[i][j];
        }
    }
    return std::sqrt(prod);
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
template<Int N, class T>
T trace(const MatXxX<N, T>& mat) {
    T prod = T(0);
    for(Int i = 0; i < N; ++i) {
        prod += mat[i][i];
    }
    return prod;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
MatXxX<N, T> dev(const MatXxX<N, T>& mat) {
    return mat - MatXxX<N, T>(LinaHelpers::trace<T>(mat) / T(N));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
void sumToDiag(MatXxX<N, T>& mat, T c) {
    for(Int i = 0; i < N; ++i) {
        mat[i][i] += c;
    }
}

template<Int N, class T>
void sum1ToDiag(MatXxX<N, T>& mat) {
    for(Int i = 0; i < N; ++i) {
        mat[i][i] += T(1.0);
    }
}

template<Int N, class T>
MatXxX<N, T> getDiagSum(const MatXxX<N, T>& mat, T c) {
    auto result = mat;
    for(Int i = 0; i < N; ++i) {
        result[i][i] += c;
    }
    return result;
}

template<Int N, class T>
MatXxX<N, T> getDiagSum(const MatXxX<N, T>& mat, const VecX<N, T>& c) {
    auto result = mat;
    for(Int i = 0; i < N; ++i) {
        result[i][i] += c[i];
    }
    return result;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
VecX<N, T> extractDiag(const MatXxX<N, T>& mat) {
    VecX<N, T> diag;
    for(Int i = 0; i < N; ++i) {
        diag[i] = mat[i][i];
    }
    return diag;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
MatXxX<N, T> diagMatrix(const VecX<N, T>& diag) {
    MatXxX<N, T> mat(0);
    for(Int i = 0; i < N; ++i) {
        mat[i][i] = diag[i];
    }
    return mat;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
void diagProduct(MatXxX<N, T>& mat, const VecX<N, T>& vec) {
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            mat[i][j] *= vec[i];
        }
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//Matrix * Matrix^-1
template<Int N, class T>
void diagProductInv(MatXxX<N, T>& mat, const VecX<N, T>& vec) {
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            mat[i][j] /= vec[i];
        }
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
T frobeniusInnerProduct(const MatXxX<N, T>& m1, const MatXxX<N, T>& m2) {
    T prod = T(0);
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            prod += m1[i][j] * m2[i][j];
        }
    }
    return prod;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
VecX<N, T> innerProduct(const VecX<N, T>& vec, const MatXxX<N, T>& mat) {
    VecX<N, T> prod(0);
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            prod[i] += vec[j] * mat[j][i];
        }
    }
    return prod;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
MatXxX<N, T> innerProduct(const MatXxX<N, T>& m1, const MatXxX<N, T>& m2) {
    MatXxX<N, T> prod(0);
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
template<class T>
Mat2x2<T> cofactor(const Mat2x2<T>& mat) {
    return Mat2x2<T>(mat[1][1], -mat[0][1],
                     -mat[1][0], mat[0][0]);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// TODO: test value row-col major
template<class T>
T elementCofactor(const Mat3x3<T>& mat, Int x, Int y) {
    T         cofactor_v;
    T         minor;
    Mat2x2<T> minor_mat;

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
template<class T>
Mat3x3<T> cofactor(const Mat3x3<T>& mat) {
    Mat2x2<T> result;
    for(Int i = 0; i < 2; ++i) {
        for(Int j = 0; j < 2; ++j) {
            result[i][j] = elementCofactor(mat, i, j);
        }
    }
    return result;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
T vonMisesPlaneStress(const Mat3x3<T>& mat) {
    const T vm = mat[0][0] * mat[0][0] + mat[1][1] * mat[1][1] + mat[2][2] * mat[2][2] -
                 mat[0][0] * mat[1][1] - mat[1][1] * mat[2][2] - mat[2][2] * mat[0][0] +
                 (mat[0][1] * mat[1][0] + mat[1][2] * mat[2][1] + mat[2][0] * mat[0][2]) * T(3.0);

    return vm > 0 ? std::sqrt(vm) : 0;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
MatXxX<N, T> randMatrix(T minVal = T(0), T maxVal = T(1.0)) {
    std::random_device                rd;
    std::mt19937                      gen(rd());
    std::uniform_real_distribution<T> dis(minVal, maxVal);

    MatXxX<N, T> result;
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            result[i][j] = dis(gen);
        }
    }
    return result;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class SizeType, class T>
StdVT<MatXxX<N, T>> randVecMatrices(SizeType size, T minVal = T(0), T maxVal = T(1.0)) {
    std::random_device                rd;
    std::mt19937                      gen(rd());
    std::uniform_real_distribution<T> dis(minVal, maxVal);

    StdVT<MatXxX<N, T>> results;
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
template<Int N, class T>
auto orientedSVD(const MatXxX<N, T>& M) {
    MatXxX<N, T> U, Vt;
    VecX<N, T>   S;

    QRSVD::svd(M, U, S, Vt);
    Vt = glm::transpose(Vt);

    MatXxX<N, T> J = MatXxX<N, T>(1.0);
    J[N - 1][N - 1] = T(-1.0);

    // Check for inversion
    if(glm::determinant(U) < T(0)) {
        U         = U * J;
        S[N - 1] *= T(-1.0);
    }
    if(glm::determinant(Vt) < T(0)) {
        Vt        = J * Vt;
        S[N - 1] *= T(-1.0);
    }

    return std::make_tuple(U, S, Vt);
} // end oriented svd

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
void QRDifferential(const Mat2x2<T>& Q, const Mat2x2<T>& R, const Mat2x2<T>& dF, Mat2x2<T>& dQ, Mat2x2<T>& dR) {
    __NT_REQUIRE(R[0][0] != 0);
    Mat2x2<T> QtdF = glm::transpose(Q) * dF;
    T         a    = -QtdF[0][1] / R[0][0];
    Mat2x2<T> QtdQ(0, -a, a, 0);
    dQ = Q * QtdQ;
    dR = glm::transpose(Q) * dF - QtdQ * R;
}

template<class T>
void QRDifferential(const Mat3x3<T>& Q, const Mat3x3<T>& R, const Mat3x3<T>& dF, Mat3x3<T>& dQ, Mat3x3<T>& dR) {
    __NT_REQUIRE(R[0][0] != 0 && R[1][1] != 0);
    Mat3x3<T> QtdF = glm::transpose(Q) * dF;
    T         w3   = QtdF[0][1] / R[0][0];
    T         w2   = -QtdF[0][2] / R[0][0];
    T         w1   = (QtdF[1][2] + w2 * R[1][0]) / R[1][1];
    Mat3x3<T> QtdQ(0, w3, -w2, -w3, 0, w1, w2, -w1, 0);
    dQ = Q * QtdQ;
    dR = glm::transpose(Q) * dF - QtdQ * R;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
auto symmetryDecomposition(const MatXxX<N, T>& M) {
    MatXxX<N, T> symComp, skewSymComp;
    auto         Mt = glm::transpose(M);

    symComp     = T(0.5) * (M + Mt);
    skewSymComp = T(0.5) * (M - Mt);
    return std::make_tuple(symComp, skewSymComp);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
auto extractFiberCotangentStress(const MatXxX<N, T>& VP, const MatXxX<N, T>& F) {
    if constexpr (N == 2) {
        return glm::outerProduct(VP[1], F[1]);
    } else {
        return MatMxN<3, 2, T>(VP[1], VP[2]) * glm::transpose(MatMxN<3, 2, T>(F[1], F[2]));
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
MatXxX<N, T> getOrthogonalSystem(const VecX<N, T>& d1) {
    if constexpr (N == 2) {
        MatXxX<N, T> M;
        M[0] = glm::normalize(d1);
        M[1] = glm::normalize(Vec2<T>(-d1.y, d1.x));
        return M;
    } else {
        const Vec3<T> basis[3] = { Vec3<T>(1, 0, 0), Vec3<T>(0, 1, 0), Vec3<T>(0, 0, 1) };
        Vec3<T>       u;
        MatXxX<N, T>  M;
        M[0] = glm::normalize(d1);

        for(UInt i = 0; i < 3; ++i) {
            if(auto de = glm::dot(M[0], basis[i]); de > T(0.1)) {
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
    return Mat4x4<typename Quaternion::value_type>(q.w, q.z, -q.y, -q.x, // col 0
                                                   -q.z, q.w, q.x, -q.y, // col 1
                                                   q.y, -q.x, q.w, -q.z, // col 2
                                                   q.x, q.y, q.z, q.w);  // col 3
}

template<class Quaternion>
auto quaternionMatrixTransposed3x4(const Quaternion& q) {
    return MatMxN<3, 4, typename Quaternion::value_type>(q.w, -q.z, q.y,    // col 0
                                                         q.z, q.w, -q.x,    // col 1
                                                         -q.y, q.x, q.w,    // col 2
                                                         -q.x, -q.y, -q.z); // col 3
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T> inline Vec4<T> quatToVec4(const Quat<T>& q) { return Vec4<T>(q.x, q.y, q.z, q.w); }
template<class T> inline Quat<T> vec4ToQuat(const Vec4<T>& q) { return Quat<T>(q.w, q.x, q.y, q.z); }

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
auto rotationMatrix2D(T angle) {
    auto      s = std::sin(angle);
    auto      c = std::cos(angle);
    Mat2x2<T> m;
    m[0][0] = c;
    m[0][1] = s;
    m[1][0] = -s;
    m[1][1] = c;
    return m;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Rotate a matrix 2 rows/N cols counter-clockwise by 90 degree
template<class T>
auto rotateRightAngle2D(const Vec2<T>& v) {
    return Vec2<T>(-v[1], v[0]);
}

template<Int N, class T>
auto rotateRightAngle2D(const glm::mat<N, 2, T>& m) {
    glm::mat<N, 2, T> result;
    for(Int i = 0; i < N; ++i) {
        result[i] = Vec2<T>(-m[i][1], m[i][0]);
    }
    return result;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// detail: https://en.wikipedia.org/wiki/Skew-symmetric_matrix#Cross_product
template<Int N, class T>
auto skewSymmetricMatrix3D(const VecX<N, T>& v) {
    return Mat3x3<T>(
        /*m[0][0] =*/ 0,
        /*m[0][1] =*/ v[2],
        /*m[0][2] =*/ -v[1],

        /*m[1][0] =*/ -v[2],
        /*m[1][1] =*/ 0,
        /*m[1][2] =*/ v[0],

        /*m[2][0] =*/ v[1],
        /*m[2][1] =*/ -v[0],
        /*m[2][2] =*/ 0
        );
}

template<Int N, class T>
auto transposedSkewSymmetricMatrix3D(const VecX<N, T>& v) {
    return Mat3x3<T>(
        /*m[0][0] =*/ 0,
        /*m[0][1] =*/ -v[2],
        /*m[0][2] =*/ v[1],

        /*m[1][0] =*/ v[2],
        /*m[1][1] =*/ 0,
        /*m[1][2] =*/ -v[0],

        /*m[2][0] =*/ -v[1],
        /*m[2][1] =*/ v[0],
        /*m[2][2] =*/ 0
        );
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
/**
 * \brief Parallel-transports u along the t1->t2 transformation.
 *
 * \param [in] u vector to be parallel-transported.
 * \param [in] t1 first vector defining transformation. It doesn't have to be normalized.
 * \param [in] t2 second vector defining transformation. It doesn't have to be normalized.
 */
template<class T>
Vec3<T> parallelTransport(const Vec3<T>& u, const Vec3<T>& t1, const Vec3<T>& t2) {
    auto b = glm::cross(t1, t2); // binormal
    if(const auto l2b = glm::length2(b); l2b < T(1e-14)) {
        return u;
    } else {
        b /= std::sqrt(l2b); // normalize b
    }
    const auto n1 = glm::normalize(glm::cross(t1, b));
    const auto n2 = glm::normalize(glm::cross(t2, b));
    return glm::dot(u, glm::normalize(t1)) * glm::normalize(t2) + glm::dot(u, n1) * n2 + glm::dot(u, b) * b;
}

/**
 * \brief Parallel-transports u along the t1->t2 transformation, assuming that u is normal to t1.
 *
 * \param [in] u vector to be parallel-transported.
 * \param [in] t1 first vector defining transformation. It doesn't have to be normalized.
 * \param [in] t2 second vector defining transformation. It doesn't have to be normalized.
 */
template<class T>
Vec3<T> normalParallelTransport(const Vec3<T>& u, const Vec3<T>& t1, const Vec3<T>& t2) {
    auto b = glm::cross(t1, t2); // binormal
    if(const auto l2b = glm::length2(b); l2b < T(1e-14)) {
        return u;
    } else {
        b /= std::sqrt(l2b); // normalize b
    }
    const auto n1 = glm::normalize(glm::cross(t1, b));
    const auto n2 = glm::normalize(glm::cross(t2, b));
    return glm::dot(u, n1) * n2 + glm::dot(u, b) * b;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
/**
 * \brief Solve system of linear equations { a1 * x + b1 * y + c1 = 0, a2 * x + b2 * y + c2 = 0 }
 */
template<class T>
std::pair<T, T> solveLinearSystem(T a1, T b1, T c1,
                                  T a2, T b2, T c2) {
    T delta = a1 * b2 - a2 * b1;
    T x     = (b1 * c2 - b2 * c1) / delta;
    T y     = (c1 * a2 - c2 * a1) / delta;
    return { x, y };
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace NTCodeBase::LinaHelpers
