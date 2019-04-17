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

#include <LibCommon/LinearAlgebra/LinearSolvers/PCGSolver.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NTCodeBase {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
void SparseColumnLowerFactor<Real_t>::clear(void) {
    nRows = 0;
    invDiag.clear();
    colValue.clear();
    colIndex.clear();
    colStart.clear();
    aDiag.clear();
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
void SparseColumnLowerFactor<Real_t>::reserve(UInt size) {
    invDiag.reserve(size);
    aDiag.reserve(size);
    colStart.reserve(size + 1);
}

template<class Real_t>
void SparseColumnLowerFactor<Real_t>::resize(UInt newSize) {
    nRows = newSize;
    invDiag.resize(nRows);
    aDiag.resize(nRows);
    colStart.resize(nRows + 1);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// PCGSolver
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
void PCGSolver<Real_t>::setSolverParameters(Real_t tolerancem_ICCPrecond, int maxIterations, Real_t MICCL0Param /*= 0.97*/, Real_t minDiagonalRatio /*= 0.25*/) {
    m_ToleranceFactor  = fmax(tolerancem_ICCPrecond, Real_t(1e-30));
    m_MaxIterations    = maxIterations;
    m_MICCL0Param      = MICCL0Param;
    m_MinDiagonalRatio = minDiagonalRatio;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
void PCGSolver<Real_t>::reserve(UInt size) {
    s.reserve(size);
    z.reserve(size);
    r.reserve(size);

    m_FixedSparseMatrix.reserve(size);
    m_ICCPrecond.reserve(size);
}

template<class Real_t>
void PCGSolver<Real_t>::resize(UInt size) {
    s.resize(size);
    z.resize(size);
    r.resize(size);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
bool PCGSolver<Real_t>::solve(const SparseMatrix<Real_t>& matrix, const StdVT<Real_t>& rhs, StdVT<Real_t>& result) {
    resize(matrix.nRows);

    result.resize(matrix.nRows);
    if(m_bZeroInitial) {
        result.assign(result.size(), 0);
    }

    m_FixedSparseMatrix.constructFromSparseMatrix(matrix);
    r = rhs;

    m_OutResidual = ParallelSTL::maxAbs<Real_t>(r);

    if(m_OutResidual < std::numeric_limits<Real_t>::min()) {
        m_OutIterations = 0;
        return true;
    }

    Real_t tol = m_ToleranceFactor * m_OutResidual;
    Real_t rho = ParallelBLAS::dotProduct<Real_t>(r, r);

    if(rho < std::numeric_limits<Real_t>::min() || std::isnan(rho)) {
        m_OutIterations = 0;
        return false;
    }

    ////////////////////////////////////////////////////////////////////////////////
    z = r;
    for(UInt iteration = 0; iteration < m_MaxIterations; ++iteration) {
        FixedSparseMatrix<Real_t>::multiply(m_FixedSparseMatrix, z, s);
        Real_t tmp = ParallelBLAS::dotProduct<Real_t>(s, z);
        if(tmp < std::numeric_limits<Real_t>::min()) {
            m_OutIterations = iteration + 1;
            return true;
        }
        Real_t alpha = rho / tmp;
        ParallelBLAS::addScaled<Real_t>(alpha,  z, result);
        ParallelBLAS::addScaled<Real_t>(-alpha, s, r);

        m_OutResidual = ParallelSTL::maxAbs<Real_t>(r);
        if(m_OutResidual < tol) {
            m_OutIterations = iteration + 1;
            return true;
        }

        Real_t rho_new = ParallelBLAS::dotProduct<Real_t>(r, r);
        Real_t beta    = rho_new / rho;
        ParallelBLAS::scaledAdd<Real_t, Real_t>(beta, r, z);
        rho = rho_new;
    }

    m_OutIterations = m_MaxIterations;
    return false;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
bool PCGSolver<Real_t>::solve_precond(const SparseMatrix<Real_t>& matrix, const StdVT<Real_t>& rhs, StdVT<Real_t>& result) {
    resize(matrix.nRows);

    result.resize(matrix.nRows);
    if(m_bZeroInitial) {
        result.assign(result.size(), 0);
    }

    m_FixedSparseMatrix.constructFromSparseMatrix(matrix);
    r = rhs;

    m_OutResidual = ParallelSTL::maxAbs<Real_t>(r);
    if(m_OutResidual < std::numeric_limits<Real_t>::min()) {
        m_OutIterations = 0;
        return true;
    }

    Real_t tol = m_ToleranceFactor * m_OutResidual;
    formPreconditioner(matrix);
    applyPreconditioner(r, z);

    Real_t rho = ParallelBLAS::dotProduct<Real_t>(z, r);
    if(rho < std::numeric_limits<Real_t>::min() || std::isnan(rho)) {
        m_OutIterations = 0;
        return false;
    }

    ////////////////////////////////////////////////////////////////////////////////
    s = z;
    for(UInt iteration = 0; iteration < m_MaxIterations; ++iteration) {
        FixedSparseMatrix<Real_t>::multiply(m_FixedSparseMatrix, s, z);
        Real_t tmp = ParallelBLAS::dotProduct<Real_t>(s, z);
        if(tmp < std::numeric_limits<Real_t>::min()) {
            m_OutIterations = iteration + 1;
            return true;
        }
        Real_t alpha = rho / tmp;
        ParallelBLAS::addScaled<Real_t>(alpha,  s, result);
        ParallelBLAS::addScaled<Real_t>(-alpha, z, r);

        m_OutResidual = ParallelSTL::maxAbs<Real_t>(r);
        if(m_OutResidual < tol) {
            m_OutIterations = iteration + 1;
            return true;
        }

        applyPreconditioner(r, z);
        Real_t rho_new = ParallelBLAS::dotProduct<Real_t>(z, r);
        Real_t beta    = rho_new / rho;
        ParallelBLAS::addScaled<Real_t, Real_t>(beta, s, z);
        s.swap(z); // s=beta*s+z
        rho = rho_new;
    }

    m_OutIterations = m_MaxIterations;
    return false;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
void PCGSolver<Real_t>::formPreconditioner(const SparseMatrix<Real_t>& matrix) {
    switch(m_PreconditionerType) {
        case Preconditioner::JACOBI:
            formPreconditioner_Jacobi(matrix);
            break;

        case Preconditioner::MICCL0:
            formPreconditioner_MICC0L0(matrix);
            break;

        case Preconditioner::MICCL0_SYMMETRIC:
            formPreconditioner_Symmetric_MICC0L0(matrix);
            break;
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
void PCGSolver<Real_t>::applyPreconditioner(const StdVT<Real_t>& x, StdVT<Real_t>& result) {
    if(m_PreconditionerType != Preconditioner::JACOBI) {
        solveLower(x, result);
        solveLower_TransposeInPlace(result);
    } else {
        applyJacobiPreconditioner(x, result);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
void PCGSolver<Real_t>::applyJacobiPreconditioner(const StdVT<Real_t>& x, StdVT<Real_t>& result) {
    ParallelExec::run<size_t>(0, x.size(),
                              [&](size_t i) {
                                  result[i] = m_JacobiPrecond[i] * x[i];
                              });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Solution routines with lower triangular matrix.
// solve L*result=rhs
template<class Real_t>
void PCGSolver<Real_t>::solveLower(const StdVT<Real_t>& rhs, StdVT<Real_t>& result) {
    result = rhs;

    for(UInt i = 0, iEnd = m_ICCPrecond.nRows; i < iEnd; ++i) {
        result[i] *= m_ICCPrecond.invDiag[i];

        for(UInt j = m_ICCPrecond.colStart[i], jEnd = m_ICCPrecond.colStart[i + 1]; j < jEnd; ++j) {
            result[m_ICCPrecond.colIndex[j]] -= m_ICCPrecond.colValue[j] * result[i];
        }
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// solve L^T*result=rhs
template<class Real_t>
void PCGSolver<Real_t>::solveLower_TransposeInPlace(StdVT<Real_t>& x) {
    UInt i = m_ICCPrecond.nRows;

    do {
        --i;
        for(UInt j = m_ICCPrecond.colStart[i], jEnd = m_ICCPrecond.colStart[i + 1]; j < jEnd; ++j) {
            x[i] -= m_ICCPrecond.colValue[j] * x[m_ICCPrecond.colIndex[j]];
        }

        x[i] *= m_ICCPrecond.invDiag[i];
    } while(i > 0);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
void PCGSolver<Real_t>::formPreconditioner_Jacobi(const SparseMatrix<Real_t>& matrix) {
    m_JacobiPrecond.resize(matrix.nRows);
    ParallelExec::run<UInt>(0, matrix.nRows,
                            [&](UInt i) {
                                UInt k = 0;
                                if(STLHelpers::Sorted::contain(matrix.colIndex[i], i, k)) {
                                    m_JacobiPrecond[i] = Real_t(1.0) / matrix.colValue[i][k];
                                } else {
                                    m_JacobiPrecond[i] = 0;
                                }
                            });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Incomplete Cholesky m_ICCPrecondization, level zero, with option for modified version.
// Set MICCL0Param between zero (regular incomplete Cholesky) and
// one (fully modified version), with values close to one usually giving the best
// results. The minDiagonalRatio parameter is used to detect and correct
// problems in m_ICCPrecondization: if a pivot is this much less than the diagonal
// entry from the original matrix, the original matrix entry is used instead.
template<class Real_t>
void PCGSolver<Real_t>::formPreconditioner_MICC0L0(const SparseMatrix<Real_t>& matrix, Real_t MICCL0Param /*= 0.97*/, Real_t minDiagonalRatio /*= 0.25*/) {
    // first copy lower triangle of matrix into m_ICCPrecond (Note: assuming A is symmetric of course!)
    m_ICCPrecond.resize(matrix.nRows);
    for(size_t i = 0, iEnd = m_ICCPrecond.invDiag.size(); i < iEnd; ++i) {
        m_ICCPrecond.invDiag[i] = 0;
        m_ICCPrecond.aDiag[i]   = 0;
    }

    m_ICCPrecond.colValue.resize(0);
    m_ICCPrecond.colIndex.resize(0);

    for(UInt i = 0, iEnd = matrix.nRows; i < iEnd; ++i) {
        m_ICCPrecond.colStart[i] = static_cast<UInt>(m_ICCPrecond.colIndex.size());

        for(UInt j = 0, jEnd = static_cast<UInt>(matrix.colIndex[i].size()); j < jEnd; ++j) {
            if(matrix.colIndex[i][j] > i) {
                m_ICCPrecond.colIndex.push_back(matrix.colIndex[i][j]);
                m_ICCPrecond.colValue.push_back(matrix.colValue[i][j]);
            } else if(matrix.colIndex[i][j] == i) {
                m_ICCPrecond.invDiag[i] = m_ICCPrecond.aDiag[i] = matrix.colValue[i][j];
            }
        }
    }

    m_ICCPrecond.colStart[matrix.nRows] = static_cast<UInt>(m_ICCPrecond.colIndex.size());

    for(UInt k = 0, kEnd = matrix.nRows; k < kEnd; ++k) {
        if(m_ICCPrecond.aDiag[k] < std::numeric_limits<Real_t>::min()) {
            continue; // null row/column
        }

        // figure out the final L(k,k) entry
        if(m_ICCPrecond.invDiag[k] < minDiagonalRatio * m_ICCPrecond.aDiag[k]) {
            m_ICCPrecond.invDiag[k] = Real_t(1.0) / std::sqrt(m_ICCPrecond.aDiag[k]); // drop to Gauss-Seidel here if the pivot looks dangerously small
        } else {
            m_ICCPrecond.invDiag[k] = Real_t(1.0) / std::sqrt(m_ICCPrecond.invDiag[k]);
        }

        // finalize the k'th column L(:,k)
        for(UInt p = m_ICCPrecond.colStart[k], pEnd = m_ICCPrecond.colStart[k + 1]; p < pEnd; ++p) {
            m_ICCPrecond.colValue[p] *= m_ICCPrecond.invDiag[k];
        }

        // incompletely eliminate L(:,k) from future columns, modifying diagonals
        for(UInt p = m_ICCPrecond.colStart[k], pEnd = m_ICCPrecond.colStart[k + 1]; p < pEnd; ++p) {
            UInt   j          = m_ICCPrecond.colIndex[p];         // work on column j
            Real_t multiplier = m_ICCPrecond.colValue[p];
            Real_t missing    = 0;
            UInt   a          = m_ICCPrecond.colStart[k];
            // first look for contributions to missing from dropped entries above the diagonal in column j
            UInt b = 0;

            while(a < m_ICCPrecond.colStart[k + 1] && m_ICCPrecond.colIndex[a] < j) {
                // look for m_ICCPrecond.rowindex[a] in matrix.index[j] starting at b
                while(b < matrix.colIndex[j].size()) {
                    if(matrix.colIndex[j][b] < m_ICCPrecond.colIndex[a]) {
                        ++b;
                    } else if(matrix.colIndex[j][b] == m_ICCPrecond.colIndex[a]) {
                        break;
                    } else {
                        missing += m_ICCPrecond.colValue[a];
                        break;
                    }
                }

                ++a;
            }

            // adjust the diagonal j,j entry
            if(a < m_ICCPrecond.colStart[k + 1] && m_ICCPrecond.colIndex[a] == j) {
                m_ICCPrecond.invDiag[j] -= multiplier * m_ICCPrecond.colValue[a];
            }

            ++a;
            // and now eliminate from the nonzero entries below the diagonal in column j (or add to missing if we can't)
            b = m_ICCPrecond.colStart[j];

            while(a < m_ICCPrecond.colStart[k + 1] && b < m_ICCPrecond.colStart[j + 1]) {
                if(m_ICCPrecond.colIndex[b] < m_ICCPrecond.colIndex[a]) {
                    ++b;
                } else if(m_ICCPrecond.colIndex[b] == m_ICCPrecond.colIndex[a]) {
                    m_ICCPrecond.colValue[b] -= multiplier * m_ICCPrecond.colValue[a];
                    ++a;
                    ++b;
                } else {
                    missing += m_ICCPrecond.colValue[a];
                    ++a;
                }
            }

            // and if there's anything left to do, add it to missing
            while(a < m_ICCPrecond.colStart[k + 1]) {
                missing += m_ICCPrecond.colValue[a];
                ++a;
            }

            // and do the final diagonal adjustment from the missing entries
            m_ICCPrecond.invDiag[j] -= MICCL0Param * multiplier * missing;
        }
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
void PCGSolver<Real_t>::formPreconditioner_Symmetric_MICC0L0(const SparseMatrix<Real_t>& matrix, Real_t minDiagonalRatio /*= 0.25*/) {
    // first copy lower triangle of matrix into m_ICCPrecond (Note: assuming A is symmetric of course!)
    m_ICCPrecond.resize(matrix.nRows);
    for(size_t i = 0, iEnd = m_ICCPrecond.invDiag.size(); i < iEnd; ++i) {
        m_ICCPrecond.invDiag[i] = 0;
        m_ICCPrecond.aDiag[i]   = 0;
    }

    m_ICCPrecond.colValue.resize(0);
    m_ICCPrecond.colIndex.resize(0);

    for(UInt i = 0, iEnd = matrix.nRows; i < iEnd; ++i) {
        m_ICCPrecond.colStart[i] = static_cast<UInt>(m_ICCPrecond.colIndex.size());

        for(UInt j = 0, jEnd = static_cast<UInt>(matrix.colIndex[i].size()); j < jEnd; ++j) {
            if(matrix.colIndex[i][j] > i) {
                m_ICCPrecond.colIndex.push_back(matrix.colIndex[i][j]);
                m_ICCPrecond.colValue.push_back(matrix.colValue[i][j]);
            } else if(matrix.colIndex[i][j] == i) {
                m_ICCPrecond.invDiag[i] = m_ICCPrecond.aDiag[i] = matrix.colValue[i][j];
            }
        }
    }

    m_ICCPrecond.colStart[matrix.nRows] = static_cast<UInt>(m_ICCPrecond.colIndex.size());
    // now do the incomplete m_ICCPrecondization (figure out numerical values)

    for(UInt k = 0, kEnd = matrix.nRows; k < kEnd; ++k) {
        if(m_ICCPrecond.aDiag[k] < std::numeric_limits<Real_t>::min()) {
            continue; // null row/column
        }

        // figure out the final L(k,k) entry
        if(m_ICCPrecond.invDiag[k] < minDiagonalRatio * m_ICCPrecond.aDiag[k]) {
            m_ICCPrecond.invDiag[k] = Real_t(1.0) / std::sqrt(m_ICCPrecond.aDiag[k]); // drop to Gauss-Seidel here if the pivot looks dangerously small
        } else {
            m_ICCPrecond.invDiag[k] = Real_t(1.0) / std::sqrt(m_ICCPrecond.invDiag[k]);
        }

        // finalize the k'th column L(:,k)
        for(UInt p = m_ICCPrecond.colStart[k], pEnd = m_ICCPrecond.colStart[k + 1]; p < pEnd; ++p) {
            m_ICCPrecond.colValue[p] *= m_ICCPrecond.invDiag[k];
        }

        for(UInt p = m_ICCPrecond.colStart[k], pEnd = m_ICCPrecond.colStart[k + 1]; p < pEnd; ++p) {
            UInt   j          = m_ICCPrecond.colIndex[p];         // work on column j
            Real_t multiplier = m_ICCPrecond.colValue[p];

            m_ICCPrecond.invDiag[j] -= multiplier * m_ICCPrecond.colValue[p];

            UInt a  = p + 1;
            UInt b  = m_ICCPrecond.colStart[j];
            UInt as = m_ICCPrecond.colStart[k + 1];
            UInt bs = m_ICCPrecond.colStart[j + 1];

            while(a < as && b < bs) {
                if(m_ICCPrecond.colIndex[b] == m_ICCPrecond.colIndex[a]) {
                    m_ICCPrecond.colValue[b] -= multiplier * m_ICCPrecond.colValue[a];
                    ++a;
                    ++b;
                } else if(m_ICCPrecond.colIndex[b] < m_ICCPrecond.colIndex[a]) {
                    ++b;
                } else {
                    ++a;
                }
            }
        }
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
__NT_INSTANTIATE_CLASS_COMMON_TYPES(PCGSolver)
__NT_INSTANTIATE_STRUCT_COMMON_TYPES(SparseColumnLowerFactor)
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace NTCodeBase
