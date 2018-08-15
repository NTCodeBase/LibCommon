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

#include <LibCommon/CommonSetup.h>
#include <LibCommon/Utils/STLHelpers.h>
#include <LibCommon/Utils/FileHelpers.h>
#include <LibCommon/Utils/NumberHelpers.h>
#include <LibCommon/ParallelHelpers/Scheduler.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Dynamic compressed sparse row matrix
//
template<class RealType>
struct SparseMatrix
{
public:
    UInt nRows;

    // for each row, a list of all column indices (sorted)
    StdVT<StdVT_UInt> colIndex;

    // values corresponding to indices
    StdVT_StdVec<RealType> colValue;

    ////////////////////////////////////////////////////////////////////////////////
    explicit SparseMatrix(UInt size = 0) : nRows(size), colIndex(size), colValue(size) {}

    void reserve(UInt size);
    void resize(UInt newSize);
    void clear();

    template<class IndexType> RealType operator()(IndexType i, IndexType j) const;
    template<class IndexType> void     setElement(IndexType i, IndexType j, RealType newValue);
    template<class IndexType> void     addElement(IndexType i, IndexType j, RealType incrementValue);
    template<class IndexType> void     eraseElement(IndexType i, IndexType j);

    void printDebug(UInt maxRows = 0) const noexcept;
    void checkSymmetry(RealType threshold = RealType(1e-8)) const noexcept;
    void printTextFile(const char* fileName);

    ////////////////////////////////////////////////////////////////////////////////
    /*static void multiply(const SparseMatrix<RealType>& matrix, const StdVT<RealType>& x, StdVT<RealType>& result);
       static void multiply_and_subtract(const SparseMatrix<RealType>& matrix, const StdVT<RealType>& x, StdVT<RealType>& result);*/
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Fixed version of SparseMatrix. This can be significantly faster for matrix-vector
// multiplies due to better data locality.
template<class RealType>
struct FixedSparseMatrix
{
    UInt nRows;

    // nonzero values row by row
    StdVT<RealType> colValue;

    // corresponding column indices
    StdVT_UInt colIndex;

    // where each row starts in value and col index (and last entry is one past the end, the number of non zeros)
    StdVT_UInt rowStart;

    ////////////////////////////////////////////////////////////////////////////////
    explicit FixedSparseMatrix(UInt size = 0) : nRows(size), colValue(0), colIndex(0), rowStart(size + 1) {}

    void reserve(UInt size) { rowStart.reserve(size + 1); }
    void resize(UInt newSize) { nRows = newSize; rowStart.resize(nRows + 1); }
    void clear(void) { colValue.resize(0); colIndex.resize(0); rowStart.resize(0); }
    void constructFromSparseMatrix(const SparseMatrix<RealType>& fixedMatrix);

    ////////////////////////////////////////////////////////////////////////////////
    static void multiply(const FixedSparseMatrix<RealType>& matrix, const StdVT<RealType>& x, StdVT<RealType>& result);
};
