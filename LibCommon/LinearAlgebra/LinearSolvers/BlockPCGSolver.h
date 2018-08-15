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

#include <LibCommon/LinearAlgebra/SparseMatrix/BlockSparseMatrix.h>
#include <LibCommon/ParallelHelpers/ParallelBLAS.h>
#include <LibCommon/ParallelHelpers/ParallelSTL.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
class BlockPCGSolver
{
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
public:
    BlockPCGSolver() = default;

    RealType residual() const noexcept { return m_OutResidual; }
    UInt     iterations() const noexcept { return m_OutIterations; }

    ////////////////////////////////////////////////////////////////////////////////
    void setSolverParameters(RealType toleranceFactor, UInt maxIterations) { m_ToleranceFactor = toleranceFactor; m_MaxIterations = maxIterations; }
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

    RealType m_ToleranceFactor = RealType(1e-20);
    UInt     m_MaxIterations   = 10000u;
    bool     m_bZeroInitial    = true;

    ////////////////////////////////////////////////////////////////////////////////
    // output
    RealType m_OutResidual   = 0;
    UInt     m_OutIterations = 0;
};
