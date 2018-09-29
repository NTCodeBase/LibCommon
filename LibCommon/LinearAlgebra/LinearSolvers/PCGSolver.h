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
template<class Real_t>
struct SparseColumnLowerFactor
{
    UInt            nRows;
    StdVT<Real_t> invDiag;  // reciprocals of diagonal elements
    StdVT_UInt      colIndex; // a list of all row indices, for each column in turn
    StdVT<Real_t> colValue; // values below the diagonal, listed column by column
    StdVT_UInt      colStart; // where each column begins in row index (plus an extra entry at the end, of #nonzeros)
    StdVT<Real_t> aDiag;    // just used in factorization: minimum "safe" diagonal entry allowed

    explicit SparseColumnLowerFactor(UInt size = 0) : nRows(size), invDiag(size), colStart(size + 1), aDiag(size) {}

    void clear(void);
    void reserve(UInt size);
    void resize(UInt newSize);
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
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
    Real_t residual() const noexcept { return m_OutResidual; }
    UInt     iterations() const noexcept { return m_OutIterations; }

    ////////////////////////////////////////////////////////////////////////////////
    void setPreconditioners(Preconditioner precond) { m_PreconditionerType = precond; }
    void setZeroInitial(bool bZeroInitial) { m_bZeroInitial = bZeroInitial; }
    void enableZeroInitial() { m_bZeroInitial = true; }
    void disableZeroInitial() { m_bZeroInitial = false; }
    void setSolverParameters(Real_t toleranceFactor, int maxIterations, Real_t MICCL0Param = Real_t(0.97), Real_t minDiagonalRatio = Real_t(0.25));
    bool solve(const SparseMatrix<Real_t>& matrix, const StdVT<Real_t>& rhs, StdVT<Real_t>& result);
    bool solve_precond(const SparseMatrix<Real_t>& matrix, const StdVT<Real_t>& rhs, StdVT<Real_t>& result);

private:
    void resize(UInt size);
    void formPreconditioner(const SparseMatrix<Real_t>& matrix);
    void applyPreconditioner(const StdVT<Real_t>& x, StdVT<Real_t>& result);
    void applyJacobiPreconditioner(const StdVT<Real_t>& x, StdVT<Real_t>& result);

    void solveLower(const StdVT<Real_t>& rhs, StdVT<Real_t>& result);
    void solveLower_TransposeInPlace(StdVT<Real_t>& x);

    void formPreconditioner_Jacobi(const SparseMatrix<Real_t>& matrix);
    void formPreconditioner_MICC0L0(const SparseMatrix<Real_t>& matrix, Real_t MICCL0Param = Real_t(0.97), Real_t minDiagonalRatio = Real_t(0.25));
    void formPreconditioner_Symmetric_MICC0L0(const SparseMatrix<Real_t>& matrix, Real_t minDiagonalRatio = Real_t(0.25));

    ////////////////////////////////////////////////////////////////////////////////
    // solver variables
    StdVT<Real_t>             z, s, r;
    FixedSparseMatrix<Real_t> m_FixedSparseMatrix;

    SparseColumnLowerFactor<Real_t> m_ICCPrecond;
    StdVT<Real_t>                   m_JacobiPrecond;

    ////////////////////////////////////////////////////////////////////////////////
    // solver parameters
    Preconditioner m_PreconditionerType = Preconditioner::MICCL0;
    Real_t       m_ToleranceFactor    = Real_t(1e-20);
    UInt           m_MaxIterations      = 10000;
    Real_t       m_MICCL0Param        = Real_t(0.97);
    Real_t       m_MinDiagonalRatio   = Real_t(0.25);
    bool           m_bZeroInitial       = true;

    ////////////////////////////////////////////////////////////////////////////////
    // output
    Real_t m_OutResidual   = 0;
    UInt     m_OutIterations = 0;
};
