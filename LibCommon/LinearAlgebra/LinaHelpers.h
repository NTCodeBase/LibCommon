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

#include <CommonSetup.h>
#include <Utils/MathHelpers.h>
#include <Utils/NumberHelpers.h>
#include <LinearAlgebra/ImplicitQRSVD.h>

#include <random>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace LinaHelpers
{
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
bool hasValidElements(const VecX<N, RealType>& vec)
{
    for(Int i = 0; i < N; ++i) {
        if(!NumberHelpers::isValidNumber(vec[i])) {
            return false;
        }
    }
    return true;
}

template<Int N, class RealType>
bool hasValidElements(const MatXxX<N, RealType>& mat)
{
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
template<Int N, class RealType>
RealType maxAbs(const MatXxX<N, RealType>& mat)
{
    RealType result = RealType(0);
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            result = MathHelpers::max(result, fabs(mat[i][j]));
        }
    }
    return result;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
RealType norm2(const MatXxX<N, RealType>& mat)
{
    RealType prod = RealType(0);
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            prod += mat[i][j] * mat[i][j];
        }
    }
    return sqrt(prod);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T, class S>
void fill(MatXxX<N, T>& mat, S x)
{
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            mat[i][j] = T(x);
        }
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
RealType trace(const MatXxX<N, RealType>& mat)
{
    RealType prod = RealType(0);
    for(Int i = 0; i < N; ++i) {
        prod += mat[i][i];
    }
    return prod;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
MatXxX<N, RealType> dev(const MatXxX<N, RealType>& mat)
{
    return mat - MatXxX<N, RealType>(LinaHelpers::trace<RealType>(mat) / RealType(N));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void sumToDiag(MatXxX<N, RealType>& mat, RealType c)
{
    for(Int i = 0; i < N; ++i) {
        mat[i][i] += c;
    }
}

template<Int N, class RealType>
void sum1ToDiag(MatXxX<N, RealType>& mat)
{
    for(Int i = 0; i < N; ++i) {
        mat[i][i] += RealType(1.0);
    }
}

template<Int N, class RealType>
MatXxX<N, RealType> getDiagSum(const MatXxX<N, RealType>& mat, RealType c)
{
    auto result = mat;
    for(Int i = 0; i < N; ++i) {
        result[i][i] += c;
    }
    return result;
}

template<Int N, class RealType>
MatXxX<N, RealType> getDiagSum(const MatXxX<N, RealType>& mat, const VecX<N, RealType>& c)
{
    auto result = mat;
    for(Int i = 0; i < N; ++i) {
        result[i][i] += c[i];
    }
    return result;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
VecX<N, RealType> extractDiag(const MatXxX<N, RealType>& mat)
{
    VecX<N, RealType> diag;
    for(Int i = 0; i < N; ++i) {
        diag[i] = mat[i][i];
    }
    return diag;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
MatXxX<N, RealType> diagMatrix(const VecX<N, RealType>& diag)
{
    MatXxX<N, RealType> mat(0);
    for(Int i = 0; i < N; ++i) {
        mat[i][i] = diag[i];
    }
    return mat;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void diagProduct(MatXxX<N, RealType>& mat, const VecX<N, RealType>& vec)
{
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            mat[i][j] *= vec[i];
        }
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//Matrix * Matrix^-1
template<Int N, class RealType>
void diagProductInv(MatXxX<N, RealType>& mat, const VecX<N, RealType>& vec)
{
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            mat[i][j] /= vec[i];
        }
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
RealType frobeniusInnerProduct(const MatXxX<N, RealType>& m1, const MatXxX<N, RealType>& m2)
{
    RealType prod = RealType(0);
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            prod += m1[i][j] * m2[i][j];
        }
    }
    return prod;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
VecX<N, RealType> innerProduct(const VecX<N, RealType>& vec, const MatXxX<N, RealType>& mat)
{
    VecX<N, RealType> prod(0);
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            prod[i] += vec[j] * mat[j][i];
        }
    }
    return prod;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
MatXxX<N, RealType> innerProduct(const MatXxX<N, RealType>& m1, const MatXxX<N, RealType>& m2)
{
    MatXxX<N, RealType> prod(0);
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
template<class RealType>
Mat2x2<RealType> cofactor(const Mat2x2<RealType>& mat)
{
    return Mat2x2<RealType>(mat[1][1], -mat[0][1],
                            -mat[1][0], mat[0][0]);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// TODO: test value row-col major
template<class RealType>
RealType elementCofactor(const Mat3x3<RealType>& mat, Int x, Int y)
{
    RealType         cofactor_v;
    RealType         minor;
    Mat2x2<RealType> minor_mat;

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
template<class RealType>
Mat3x3<RealType> cofactor(const Mat3x3<RealType>& mat)
{
    Mat2x2<RealType> result;
    for(Int i = 0; i < 2; ++i) {
        for(Int j = 0; j < 2; ++j) {
            result[i][j] = elementCofactor(mat, i, j);
        }
    }
    return result;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class RealType>
RealType vonMisesPlaneStress(const Mat3x3<RealType>& mat)
{
    const RealType vm = mat[0][0] * mat[0][0] + mat[1][1] * mat[1][1] + mat[2][2] * mat[2][2] -
                        mat[0][0] * mat[1][1] - mat[1][1] * mat[2][2] - mat[2][2] * mat[0][0] +
                        (mat[0][1] * mat[1][0] + mat[1][2] * mat[2][1] + mat[2][0] * mat[0][2]) * RealType(3.0);

    return vm > 0 ? std::sqrt(vm) : 0;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
MatXxX<N, RealType> randMatrix(RealType minVal = RealType(0), RealType maxVal = RealType(1.0))
{
    std::random_device                       rd;
    std::mt19937                             gen(rd());
    std::uniform_real_distribution<RealType> dis(minVal, maxVal);

    MatXxX<N, RealType> result;
    for(Int i = 0; i < N; ++i) {
        for(Int j = 0; j < N; ++j) {
            result[i][j] = dis(gen);
        }
    }
    return result;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class SizeType, class RealType>
StdVT<MatXxX<N, RealType>> randVecMatrices(SizeType size, RealType minVal = RealType(0), RealType maxVal = RealType(1.0))
{
    std::random_device                       rd;
    std::mt19937                             gen(rd());
    std::uniform_real_distribution<RealType> dis(minVal, maxVal);

    StdVT<MatXxX<N, RealType>> results;
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
template<Int N, class RealType>
auto orientedSVD(const MatXxX<N, RealType>& M)
{
    MatXxX<N, RealType> U, Vt;
    VecX<N, RealType>   S;

    QRSVD::svd(M, U, S, Vt);
    Vt = glm::transpose(Vt);

    MatXxX<N, RealType> J = MatXxX<N, RealType>(1.0);
    J[N - 1][N - 1] = RealType(-1.0);

    // Check for inversion
    if(glm::determinant(U) < RealType(0)) {
        U         = U * J;
        S[N - 1] *= RealType(-1.0);
    }
    if(glm::determinant(Vt) < RealType(0)) {
        Vt        = J * Vt;
        S[N - 1] *= RealType(-1.0);
    }

    return std::make_tuple(U, S, Vt);
}       // end oriented svd

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class RealType>
void QRDifferential(const Mat2x2<RealType>& Q, const Mat2x2<RealType>& R, const Mat2x2<RealType>& dF, Mat2x2<RealType>& dQ, Mat2x2<RealType>& dR)
{
    __NT_REQUIRE(R[0][0] != 0);
    Mat2x2<RealType> QtdF = glm::transpose(Q) * dF;
    RealType         a    = -QtdF[0][1] / R[0][0];
    Mat2x2<RealType> QtdQ(0, -a, a, 0);
    dQ = Q * QtdQ;
    dR = glm::transpose(Q) * dF - QtdQ * R;
}

template<class RealType>
void QRDifferential(const Mat3x3<RealType>& Q, const Mat3x3<RealType>& R, const Mat3x3<RealType>& dF, Mat3x3<RealType>& dQ, Mat3x3<RealType>& dR)
{
    __NT_REQUIRE(R[0][0] != 0 && R[1][1] != 0);
    Mat3x3<RealType> QtdF = glm::transpose(Q) * dF;
    RealType         w3   = QtdF[0][1] / R[0][0];
    RealType         w2   = -QtdF[0][2] / R[0][0];
    RealType         w1   = (QtdF[1][2] + w2 * R[1][0]) / R[1][1];
    Mat3x3<RealType> QtdQ(0, w3, -w2, -w3, 0, w1, w2, -w1, 0);
    dQ = Q * QtdQ;
    dR = glm::transpose(Q) * dF - QtdQ * R;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
auto symmetryDecomposition(const MatXxX<N, RealType>& M)
{
    MatXxX<N, RealType> symComp, skewSymComp;
    auto                Mt = glm::transpose(M);

    symComp     = RealType(0.5) * (M + Mt);
    skewSymComp = RealType(0.5) * (M - Mt);
    return std::make_tuple(symComp, skewSymComp);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
auto extractFiberCotangentStress(const MatXxX<N, RealType>& VP, const MatXxX<N, RealType>& F)
{
    if constexpr (N == 2) {
        return glm::outerProduct(VP[1], F[1]);
    } else {
        return MatMxN<3, 2, RealType>(VP[1], VP[2]) * glm::transpose(MatMxN<3, 2, RealType>(F[1], F[2]));
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
}   // end namespace LinaHelpers
