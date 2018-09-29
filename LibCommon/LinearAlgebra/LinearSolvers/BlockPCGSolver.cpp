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

#include <LibCommon/LinearAlgebra/LinearSolvers/BlockPCGSolver.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
bool BlockPCGSolver<N, Real_t>::solve(const BlockSparseMatrix<MatNxN>& matrix, const StdVT_VecN& rhs, StdVT_VecN& result)
{
    const UInt n = matrix.size();
    if(z.size() != n) {
        s.resize(n);
        z.resize(n);
        r.resize(n);
    }

    if(m_bZeroInitial) {
        result.assign(result.size(), VecN(0));
    }

    m_FixedSparseMatrix.constructFromSparseMatrix(matrix);
    r = rhs;

    m_OutResidual = ParallelSTL::maxAbs<N, Real_t>(r);

    if(m_OutResidual < m_ToleranceFactor) {
        m_OutIterations = 0;
        return true;
    }

    Real_t tol = m_ToleranceFactor * m_OutResidual;
    Real_t rho = ParallelBLAS::dotProduct<N, Real_t>(r, r);

    if(std::fabs(rho) < m_ToleranceFactor || std::isnan(rho)) {
        m_OutIterations = 0;
        return true;
    }

    z = r;

    for(UInt iteration = 0; iteration < m_MaxIterations; ++iteration) {
        FixedBlockSparseMatrix<MatNxN>::multiply(m_FixedSparseMatrix, z, s);
        Real_t tmp = ParallelBLAS::dotProduct<N, Real_t>(s, z);

        if(std::fabs(tmp) < m_ToleranceFactor || std::isnan(tmp)) {
            m_OutIterations = iteration;
            return true;
        }

        Real_t alpha = rho / tmp;

        tbb::parallel_invoke(
            [&] { ParallelBLAS::addScaled(alpha, z, result); },
            [&] { ParallelBLAS::addScaled(-alpha, s, r); }
            );

        m_OutResidual = ParallelSTL::maxAbs<N, Real_t>(r);

        if(m_OutResidual < tol) {
            m_OutIterations = iteration + 1;
            return true;
        }

        Real_t rho_new = ParallelBLAS::dotProduct<N, Real_t>(r, r);
        Real_t beta    = rho_new / rho;
        ParallelBLAS::scaledAdd(beta, r, z);
        rho = rho_new;
    }

    m_OutIterations = m_MaxIterations;
    return false;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
bool BlockPCGSolver<N, Real_t>::solve_precond(const BlockSparseMatrix<MatNxN>& matrix, const StdVT_VecN& rhs, StdVT_VecN& result)
{
    UInt n = matrix.size();

    if(z.size() != n) {
        s.resize(n);
        z.resize(n);
        r.resize(n);
    }

    // zero out the result
    if(m_bZeroInitial) {
        result.assign(result.size(), VecN(0));
    }

    m_FixedSparseMatrix.constructFromSparseMatrix(matrix);
    r = rhs;

    FixedBlockSparseMatrix<MatNxN>::multiply(m_FixedSparseMatrix, result, s);
    ParallelBLAS::addScaled(Real_t(-1.0), s, r);

    m_OutResidual = ParallelSTL::maxAbs<N, Real_t>(r);
    if(m_OutResidual < m_ToleranceFactor) {
        m_OutIterations = 0;
        return true;
    }

    Real_t tol = m_ToleranceFactor * m_OutResidual;

    formPreconditioner(matrix);
    applyPreconditioner(r, z);

    Real_t rho = ParallelBLAS::dotProduct<N, Real_t>(z, r);
    if(std::fabs(rho) < m_ToleranceFactor || std::isnan(rho)) {
        m_OutIterations = 0;
        return true;
    }

    s = z;

    for(UInt iteration = 0; iteration < m_MaxIterations; ++iteration) {
        FixedBlockSparseMatrix<MatNxN>::multiply(m_FixedSparseMatrix, s, z);
        Real_t tmp = ParallelBLAS::dotProduct<N, Real_t>(s, z);

        if(std::fabs(tmp) < m_ToleranceFactor || std::isnan(tmp)) {
            m_OutIterations = iteration;
            return true;
        }

        Real_t alpha = rho / tmp;

        tbb::parallel_invoke(
            [&] { ParallelBLAS::addScaled(alpha, s, result); },
            [&] { ParallelBLAS::addScaled(-alpha, z, r); }
            );

        m_OutResidual = ParallelSTL::maxAbs<N, Real_t>(r);
        if(m_OutResidual < tol) {
            m_OutIterations = iteration + 1;
            return true;
        }

        applyPreconditioner(r, z);

        Real_t rho_new = ParallelBLAS::dotProduct<N, Real_t>(z, r);
        Real_t beta    = rho_new / rho;

        ParallelBLAS::addScaled(beta, s, z);
        s.swap(z);                     // s=beta*s+z
        rho = rho_new;
    }

    m_OutIterations = m_MaxIterations;
    return false;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void BlockPCGSolver<N, Real_t>::formPreconditioner(const BlockSparseMatrix<MatNxN>& matrix)
{
    m_JacobiPreconditioner.resize(matrix.size());
    Scheduler::parallel_for(matrix.size(),
                            [&](UInt i)
                            {
                                const auto& v = matrix.getIndices(i);
                                const auto it = std::lower_bound(v.begin(), v.end(), i);
                                MatNxN tmp_inv;
                                if(it != v.end()) {
                                    MatNxN tmp = matrix.getValues(i)[std::distance(v.begin(), it)];
                                    for(Int j = 0; j < MatNxN::length(); ++j) {
                                        for(Int k = 0; k < MatNxN::length(); ++k) {
                                            tmp_inv[k][j] = (j == k && tmp[j][j] != 0) ? typename MatNxN::value_type(1.0) / tmp[j][j] : 0;
                                        }
                                    }
                                }
                                m_JacobiPreconditioner[i] = tmp_inv;
                            });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void BlockPCGSolver<N, Real_t>::applyPreconditioner(const StdVT_VecN& x, StdVT_VecN& result)
{
    Scheduler::parallel_for(x.size(), [&](size_t i) { result[i] = m_JacobiPreconditioner[i] * x[i]; });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(BlockPCGSolver)
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
