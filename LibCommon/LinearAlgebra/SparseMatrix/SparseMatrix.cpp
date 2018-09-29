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

#include <LibCommon/LinearAlgebra/SparseMatrix/SparseMatrix.h>
#include <LibCommon/Utils/Formatters.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Dynamic compressed sparse row matrix.
//
template<class Real_t>
void SparseMatrix<Real_t>::reserve(UInt size)
{
    colIndex.reserve(size);
    colValue.reserve(size);
}

template<class Real_t>
void SparseMatrix<Real_t>::resize(UInt newSize)
{
    nRows = newSize;
    colIndex.resize(nRows);
    colValue.resize(nRows);
}

template<class Real_t>
void SparseMatrix<Real_t>::clear()
{
    for(UInt i = 0; i < nRows; ++i) {
        colIndex[i].resize(0);
        colValue[i].resize(0);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
template<class IndexType>
Real_t SparseMatrix<Real_t>::operator()(IndexType i, IndexType j) const
{
    assert(static_cast<UInt>(i) < nRows && static_cast<UInt>(j) < nRows);
    UInt k = 0;
    if(STLHelpers::Sorted::contain(colIndex[i], static_cast<UInt>(j), k)) {
        return colValue[i][k];
    } else {
        return 0;
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
template<class IndexType>
void SparseMatrix<Real_t>::setElement(IndexType i, IndexType j, Real_t newValue)
{
    assert(static_cast<UInt>(i) < nRows && static_cast<UInt>(j) < nRows);
    UInt k = 0;
    if(STLHelpers::Sorted::contain(colIndex[i], static_cast<UInt>(j), k)) {
        colValue[i][k] = newValue;
    } else {
        STLHelpers::Sorted::insertPairSorted(colIndex[i], static_cast<UInt>(j), colValue[i], newValue);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
template<class IndexType>
void SparseMatrix<Real_t>::addElement(IndexType i, IndexType j, Real_t incrementValue)
{
    assert(static_cast<UInt>(i) < nRows && static_cast<UInt>(j) < nRows);
    UInt k = 0;
    if(STLHelpers::Sorted::contain(colIndex[i], static_cast<UInt>(j), k)) {
        colValue[i][k] += incrementValue;
    } else {
        STLHelpers::Sorted::insertPairSorted(colIndex[i], static_cast<UInt>(j), colValue[i], incrementValue);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
template<class IndexType>
void SparseMatrix<Real_t>::eraseElement(IndexType i, IndexType j)
{
    assert(static_cast<UInt>(i) < nRows && static_cast<UInt>(j) < nRows);
    UInt k = 0;
    if(STLHelpers::Sorted::contain(colIndex[i], static_cast<UInt>(j), k)) {
        colIndex[i].erase(colIndex[i].begin() + k);
        colValue[i].erase(colValue[i].begin() + k);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
void SparseMatrix<Real_t>::printDebug(UInt maxRows /*= 0*/) const noexcept
{
    if(maxRows == 0) {
        maxRows = nRows;
    }

    for(UInt i = 0; i < maxRows; ++i) {
        if(colIndex[i].size() == 0) {
            continue;
        }

        std::cout << "Line " << i << ": ";

        for(UInt j = 0; j < colIndex[i].size(); ++j) {
            std::cout << colIndex[i][j] << "(" << Formatters::toSciString(colValue[i][j]) << "), ";
        }

        std::cout << std::endl;
    }

    std::cout << std::endl;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
void SparseMatrix<Real_t>::checkSymmetry(Real_t threshold /* = Real_t(1e-8) */) const noexcept
{
    bool check = true;
    std::cout << "============================== Checking Matrix Symmetry... ==============================" << std::endl;
    std::cout << "Matrix size: " << nRows << std::endl;

    Scheduler::parallel_for(nRows,
                            [&](UInt i)
                            {
                                for(UInt j = i + 1; j < nRows; ++j) {
                                    if(STLHelpers::Sorted::contain(colIndex[i], j)) {
                                        auto err = std::abs((*this)(i, j) - (*this)(j, i));
                                        if(err > threshold) {
                                            check = false;
                                            std::cout << "Invalid matrix element at index " << i << ", " << j
                                                      << ", err = " << err << ": "
                                                      << "matrix(" << i << ", " << j << ") = " << (*this)(i, j) << " != "
                                                      << "matrix(" << j << ", " << i << ") = " << (*this)(j, i) << std::endl;
                                        }
                                    }
                                }
                            });

    if(check) {
        std::cout << "All matrix elements are valid!" << std::endl;
    } else {
        std::cout << "There are some invalid matrix elements!" << std::endl;
    }

    std::cout << std::endl;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
void SparseMatrix<Real_t>::printTextFile(const char* fileName)
{
    StdVT<String> matContent;

    for(UInt i = 0; i < nRows; ++i) {
        if(colIndex[i].size() == 0) {
            continue;
        }
        for(size_t j = 0; j < colIndex[i].size(); ++j) {
            String str = std::to_string(i + 1);
            str += "    ";
            str += std::to_string(colIndex[i][j] + 1);
            str += "    ";
            str += std::to_string(colValue[i][j]);
            matContent.push_back(str);
        }
    }

    FileHelpers::writeFile(matContent, fileName);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Fixed version of SparseMatrix
//
template<class Real_t>
void FixedSparseMatrix<Real_t>::constructFromSparseMatrix(const SparseMatrix<Real_t>& matrix)
{
    resize(matrix.nRows);
    rowStart[0] = 0;
    for(UInt i = 0; i < nRows; ++i) {
        rowStart[i + 1] = rowStart[i] + static_cast<UInt>(matrix.colIndex[i].size());
    }

    // in cases the matrix has empty row, accessing start row index may be out of range
    // so, add one extra element for padding
    colIndex.resize(rowStart[nRows] + 1);
    colValue.resize(rowStart[nRows] + 1);

    Scheduler::parallel_for(matrix.nRows,
                            [&](UInt i)
                            {
                                memcpy(&colIndex[rowStart[i]], matrix.colIndex[i].data(), matrix.colIndex[i].size() * sizeof(UInt));
                                memcpy(&colValue[rowStart[i]], matrix.colValue[i].data(), matrix.colValue[i].size() * sizeof(Real_t));
                            });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// perform result=matrix*x
template<class Real_t>
void FixedSparseMatrix<Real_t>::multiply(const FixedSparseMatrix<Real_t>& matrix, const StdVT<Real_t>& x, StdVT<Real_t>& result)
{
    assert(matrix.nRows == static_cast<UInt>(x.size()));
    result.resize(matrix.nRows);
    Scheduler::parallel_for(matrix.nRows,
                            [&](UInt i)
                            {
                                Real_t tmpResult = 0;
                                for(UInt j = matrix.rowStart[i], jEnd = matrix.rowStart[i + 1]; j < jEnd; ++j) {
                                    tmpResult += matrix.colValue[j] * x[matrix.colIndex[j]];
                                }
                                result[i] = tmpResult;
                            });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define __BNN_INSTANTIATE_SPARSE_MATRIX_FUNCS(IntType, Real_t)                                              \
    template Real_t SparseMatrix<Real_t >::operator()(IntType i, IntType j) const;                        \
    template void SparseMatrix<Real_t>::setElement<IntType>(IntType i, IntType j, Real_t newValue);       \
    template void SparseMatrix<Real_t>::addElement<IntType>(IntType i, IntType j, Real_t incrementValue); \
    template void SparseMatrix<Real_t>::eraseElement<IntType>(IntType i, IntType j);

__BNN_INSTANTIATE_SPARSE_MATRIX_FUNCS(Int,    float)
__BNN_INSTANTIATE_SPARSE_MATRIX_FUNCS(UInt,   float)
__BNN_INSTANTIATE_SPARSE_MATRIX_FUNCS(Int64,  float)
__BNN_INSTANTIATE_SPARSE_MATRIX_FUNCS(UInt64, float)

__BNN_INSTANTIATE_SPARSE_MATRIX_FUNCS(Int,    double)
__BNN_INSTANTIATE_SPARSE_MATRIX_FUNCS(UInt,   double)
__BNN_INSTANTIATE_SPARSE_MATRIX_FUNCS(Int64,  double)
__BNN_INSTANTIATE_SPARSE_MATRIX_FUNCS(UInt64, double)

__NT_INSTANTIATE_STRUCT_COMMON_TYPES(SparseMatrix)
__NT_INSTANTIATE_STRUCT_COMMON_TYPES(FixedSparseMatrix)
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
