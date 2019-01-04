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

#include <LibCommon/LinearAlgebra/SparseMatrix/BlockSparseMatrix.h>
#include <LibCommon/ParallelHelpers/ParallelBLAS.h>
#include <LibCommon/ParallelHelpers/ParallelSTL.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NTCodeBase {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class BlockPCGSolver {
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIAS
    ////////////////////////////////////////////////////////////////////////////////
public:
    BlockPCGSolver() = default;

    Real_t residual() const noexcept { return m_OutResidual; }
    UInt     iterations() const noexcept { return m_OutIterations; }

    ////////////////////////////////////////////////////////////////////////////////
    void setSolverParameters(Real_t toleranceFactor, UInt maxIterations) { m_ToleranceFactor = toleranceFactor; m_MaxIterations = maxIterations; }
    void setZeroInitial(bool bZeroInitial) { m_bZeroInitial = bZeroInitial; }
    void enableZeroInitial() { m_bZeroInitial = true; }
    void disableZeroInitial() { m_bZeroInitial = false; }

    bool solve(const BlockSparseMatrix<MatNxN>& matrix, const StdVT_VecN& rhs, StdVT_VecN& result);
    bool solve_precond(const BlockSparseMatrix<MatNxN>& matrix, const StdVT_VecN& rhs, StdVT_VecN& result);

private:
    void formPreconditioner(const BlockSparseMatrix<MatNxN>& matrix);
    void applyPreconditioner(const StdVT_VecN& x, StdVT_VecN& result);

    ////////////////////////////////////////////////////////////////////////////////
    StdVT_VecN                     z, s, r;
    StdVT_MatNxN                   m_JacobiPreconditioner;
    FixedBlockSparseMatrix<MatNxN> m_FixedSparseMatrix;

    Real_t m_ToleranceFactor = Real_t(1e-20);
    UInt   m_MaxIterations   = 10000u;
    bool   m_bZeroInitial    = true;

    ////////////////////////////////////////////////////////////////////////////////
    // output
    Real_t m_OutResidual   = 0;
    UInt   m_OutIterations = 0;
};
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace NTCodeBase
