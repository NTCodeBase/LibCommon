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

#include <LibCommon/LinearAlgebra/SparseMatrix/SparseMatrix.h>
#include <LibCommon/ParallelHelpers/ParallelBLAS.h>
#include <LibCommon/ParallelHelpers/ParallelSTL.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// A simple compressed sparse column data structure (with separate diagonal)
// for lower triangular matrices
template<class RealType>
struct SparseColumnLowerFactor
{
    UInt            nRows;
    StdVT<RealType> invDiag;  // reciprocals of diagonal elements
    StdVT_UInt      colIndex; // a list of all row indices, for each column in turn
    StdVT<RealType> colValue; // values below the diagonal, listed column by column
    StdVT_UInt      colStart; // where each column begins in row index (plus an extra entry at the end, of #nonzeros)
    StdVT<RealType> aDiag;    // just used in factorization: minimum "safe" diagonal entry allowed

    explicit SparseColumnLowerFactor(UInt size = 0) : nRows(size), invDiag(size), colStart(size + 1), aDiag(size) {}

    void clear(void);
    void reserve(UInt size);
    void resize(UInt newSize);
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class RealType>
class PCGSolver
{
public:
    enum Preconditioner
    {
        JACOBI,
        MICCL0,
        MICCL0_SYMMETRIC
    };

    PCGSolver() = default;
    void reserve(UInt size);
    RealType residual() const noexcept { return m_OutResidual; }
    UInt     iterations() const noexcept { return m_OutIterations; }

    ////////////////////////////////////////////////////////////////////////////////
    void setPreconditioners(Preconditioner precond) { m_PreconditionerType = precond; }
    void setZeroInitial(bool bZeroInitial) { m_bZeroInitial = bZeroInitial; }
    void enableZeroInitial() { m_bZeroInitial = true; }
    void disableZeroInitial() { m_bZeroInitial = false; }
    void setSolverParameters(RealType toleranceFactor, int maxIterations, RealType MICCL0Param = RealType(0.97), RealType minDiagonalRatio = RealType(0.25));
    bool solve(const SparseMatrix<RealType>& matrix, const StdVT<RealType>& rhs, StdVT<RealType>& result);
    bool solve_precond(const SparseMatrix<RealType>& matrix, const StdVT<RealType>& rhs, StdVT<RealType>& result);

private:
    void resize(UInt size);
    void formPreconditioner(const SparseMatrix<RealType>& matrix);
    void applyPreconditioner(const StdVT<RealType>& x, StdVT<RealType>& result);
    void applyJacobiPreconditioner(const StdVT<RealType>& x, StdVT<RealType>& result);

    void solveLower(const StdVT<RealType>& rhs, StdVT<RealType>& result);
    void solveLower_TransposeInPlace(StdVT<RealType>& x);

    void formPreconditioner_Jacobi(const SparseMatrix<RealType>& matrix);
    void formPreconditioner_MICC0L0(const SparseMatrix<RealType>& matrix, RealType MICCL0Param = RealType(0.97), RealType minDiagonalRatio = RealType(0.25));
    void formPreconditioner_Symmetric_MICC0L0(const SparseMatrix<RealType>& matrix, RealType minDiagonalRatio = RealType(0.25));

    ////////////////////////////////////////////////////////////////////////////////
    // solver variables
    StdVT<RealType>             z, s, r;
    FixedSparseMatrix<RealType> m_FixedSparseMatrix;

    SparseColumnLowerFactor<RealType> m_ICCPrecond;
    StdVT<RealType>                   m_JacobiPrecond;

    ////////////////////////////////////////////////////////////////////////////////
    // solver parameters
    Preconditioner m_PreconditionerType = Preconditioner::MICCL0;
    RealType       m_ToleranceFactor    = RealType(1e-20);
    UInt           m_MaxIterations      = 10000;
    RealType       m_MICCL0Param        = RealType(0.97);
    RealType       m_MinDiagonalRatio   = RealType(0.25);
    bool           m_bZeroInitial       = true;

    ////////////////////////////////////////////////////////////////////////////////
    // output
    RealType m_OutResidual   = 0;
    UInt     m_OutIterations = 0;
};
