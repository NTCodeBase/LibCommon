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

#include <LibCommon/CommonSetup.h>
#include <LibCommon/Array/Array.h>
#include <LibCommon/Math/MathHelpers.h>
#include <LibCommon/ParallelHelpers/ParallelObjects.h>
#include <cassert>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NTCodeBase {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class Grid {
    ////////////////////////////////////////////////////////////////////////////////
    NT_TYPE_ALIAS
    ////////////////////////////////////////////////////////////////////////////////
public:
    Grid() = default;
    Grid(const VecN& bMin, const VecN& bMax, Real_t cellSize) : m_BMin(bMin), m_BMax(bMax) { setCellSize(cellSize); }

    ////////////////////////////////////////////////////////////////////////////////
    // Setters
    void setGrid(const VecN& bMin, const VecN& bMax, Real_t cellSize, Real_t clampEdge = 0);
    void setCellSize(Real_t cellSize);

    ////////////////////////////////////////////////////////////////////////////////
    // Getters
    const auto& getBMin() const noexcept { return m_BMin; }
    const auto& getBMax() const noexcept { return m_BMax; }

    const auto& getNCells() const noexcept { return m_NCells; }
    const auto& getNNodes() const noexcept { return m_NNodes; }
    auto        getNTotalCells() const noexcept { return m_NTotalCells; }
    auto        getNTotalNodes() const noexcept { return m_NTotalNodes; }

    ////////////////////////////////////////////////////////////////////////////////
    auto getCellSize() const noexcept { return m_CellSize; }
    auto getInvCellSize() const noexcept { return m_InvCellSize; }
    auto getHalfCellSize() const noexcept { return m_HalfCellSize; }
    auto getCellSizeSquared() const noexcept { return m_CellSizeSqr; }
    auto getInvCellSizeSquared() const noexcept { return m_InvCellSizeSqr; }

    ////////////////////////////////////////////////////////////////////////////////
    auto getCellVolume() const noexcept { return m_CellVolume; }

    ////////////////////////////////////////////////////////////////////////////////
    // Grid 2D =>
    template<class IndexType>
    IndexType getFlatIndex(IndexType i, IndexType j) const {
        static_assert(N == 2, "Array dimension != 2");
        return j * static_cast<IndexType>(getNCells()[0]) + i;
    }

    template<class IndexType>
    IndexType getNodeLinearizedIndex(IndexType i, IndexType j) const {
        static_assert(N == 2, "Array dimension != 2");
        return j * static_cast<IndexType>(getNNodes()[0]) + i;
    }

    template<class IndexType>
    bool isValidCell(IndexType i, IndexType j) const noexcept{
        static_assert(N == 2, "Array dimension != 2");
        return (i >= 0 &&
                j >= 0 &&
                static_cast<UInt>(i) < m_NCells[0] &&
                static_cast<UInt>(j) < m_NCells[1]);
    }

    template<class IndexType>
    bool isValidCell(const Vec2<IndexType>& index) const noexcept{
        return isValidCell(index[0], index[1]);
    }

    template<class IndexType>
    bool isValidNode(IndexType i, IndexType j) const noexcept{
        static_assert(N == 2, "Array dimension != 2");
        return (i >= 0 &&
                j >= 0 &&
                static_cast<UInt>(i) < m_NNodes[0] &&
                static_cast<UInt>(j) < m_NNodes[1]);
    }

    template<class IndexType>
    bool isValidNode(const Vec2<IndexType>& index) const noexcept{
        return isValidNode(index[0], index[1]);
    }

    // => Grid 2D
    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////
    // Grid 3D =>
    template<class IndexType>
    IndexType getFlatIndex(IndexType i, IndexType j, IndexType k) const {
        static_assert(N == 3, "Array dimension != 3");
        return (k * static_cast<IndexType>(getNCells()[1]) + j) * static_cast<IndexType>(getNCells()[0]) + i;
    }

    template<class IndexType>
    IndexType getNodeLinearizedIndex(IndexType i, IndexType j, IndexType k) const {
        static_assert(N == 3, "Array dimension != 3");
        return (k * static_cast<IndexType>(getNNodes()[1]) + j) * static_cast<IndexType>(getNNodes()[0]) + i;
    }

    template<class IndexType>
    bool isValidCell(IndexType i, IndexType j, IndexType k) const noexcept{
        static_assert(N == 3, "Array dimension != 3");
        return (i >= 0 &&
                j >= 0 &&
                k >= 0 &&
                static_cast<UInt>(i) < m_NCells[0] &&
                static_cast<UInt>(j) < m_NCells[1] &&
                static_cast<UInt>(k) < m_NCells[2]);
    }

    template<class IndexType>
    bool isValidCell(const Vec3<IndexType>& index) const noexcept{
        return isValidCell(index[0], index[1], index[2]);
    }

    template<class IndexType>
    bool isValidNode(IndexType i, IndexType j, IndexType k) const noexcept{
        static_assert(N == 3, "Array dimension != 3");
        return (i >= 0 &&
                j >= 0 &&
                k >= 0 &&
                static_cast<UInt>(i) < m_NNodes[0] &&
                static_cast<UInt>(j) < m_NNodes[1] &&
                static_cast<UInt>(k) < m_NNodes[2]);
    }

    template<class IndexType>
    bool isValidNode(const Vec3<IndexType>& index) const noexcept{
        return isValidNode(index[0], index[1], index[2]);
    }

    // => Grid 3D
    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////
    template<class IndexType>
    auto getCellIdx(const VecN& ppos) const noexcept{
        VecX<N, IndexType> cellIdx;
        for(Int d = 0; d < N; ++d) {
            cellIdx[d] = static_cast<IndexType>((ppos[d] - m_BMin[d]) / m_CellSize);
        }
        return cellIdx;
    }

    template<class IndexType>
    auto getNearestValidCellIdx(const VecX<N, IndexType>& cellIdx) const noexcept{
        VecX<N, IndexType> validCellIdx;
        for(Int d = 0; d < N; ++d) {
            validCellIdx[d] = MathHelpers::clamp<IndexType>(cellIdx[d], 0, static_cast<IndexType>(m_NCells[d]) - 1);
        }
        return validCellIdx;
    }

    template<class IndexType>
    auto getValidCellIdx(const VecN& ppos) const noexcept{
        return getNearestValidCellIdx<IndexType>(getCellIdx<IndexType>(ppos));
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Particle processing
    template<class IndexType>
    auto getWorldCoordinate(const VecX<N, IndexType>& nodeIdx) const {
        return VecN(nodeIdx) * m_CellSize + m_BMin;
    }

    template<class IndexType>
    auto getWorldCoordinate(IndexType i, IndexType j) const {
        //        static_assert(N == 2, "Array dimension != 2");
        if constexpr (N == 2)
        return Vec2<Real_t>(i, j) * m_CellSize + m_BMin;
        else {
            return Vec2<Real_t>(0);
        }
    }

    template<class IndexType>
    auto getWorldCoordinate(IndexType i, IndexType j, IndexType k) const {
        static_assert(N == 3, "Array dimension != 3");
        return Vec3<Real_t>(i, j, k) * m_CellSize + m_BMin;
    }

    auto getGridCoordinate(const VecN& ppos) const { return (ppos - m_BMin) * m_InvCellSize; }
    void getGridCoordinate(const StdVT_VecN& positions, StdVT_VecN& gridCoordinates) const;
    ////////////////////////////////////////////////////////////////////////////////
    bool isInsideGrid(const VecN& ppos) const noexcept;
    bool isInsideClampedBoundary(const VecN& ppos) const noexcept;
    VecN constrainedBoundaryPosition(const VecN& positions) const noexcept;
    VecN constrainedClampedBoundaryPosition(const VecN& positions) const noexcept;
    void constrainedBoundaryInPlace(VecN& positions) const noexcept;
    void constrainedClampedBoundaryInPlace(VecN& positions) const noexcept;

    void constrainToGridBoundary(StdVT_VecN& positions);
    void constrainToClampedBoundary(StdVT_VecN& positions);
    void collectIndexToCells(const StdVT_VecN& positions);
    void collectIndexToCells(const StdVT_VecN& positions, StdVT<VecX<N, Int>>& particleCellIdx);
    void collectIndexToCells(const StdVT_VecN& positions, StdVT_VecN& gridCoordinates);
    void getNeighborList(const StdVT_VecN& positions, StdVT<StdVT_UInt>& neighborList, Int cellSpan = 1);
    void getNeighborList(const VecN& ppos, StdVT_UInt& neighborList, Int cellSpan = 1);
    void getNeighborList(const StdVT_VecN& positions, StdVT<StdVT_UInt>& neighborList, Real_t d2, Int cellSpan = 1);
    void getNeighborList(const StdVT_VecN& positions, const VecN& ppos, StdVT_UInt& neighborList, Real_t d2, Int cellSpan = 1);
    void sortData(StdVT_VecN& data);

    const StdVT_UInt& getParticleIdxSortedByCell();

    template<class IndexType>
    const auto& getParticleIdxInCell(const VecX<N, IndexType>& cellIdx) const { return m_ParticleIdxInCell(cellIdx); }
protected:
    VecN   m_BMin           = VecN(-1.0);
    VecN   m_BMax           = VecN(1.0);
    VecN   m_ClampedBMin    = VecN(-1.0);
    VecN   m_ClampedBMax    = VecN(1.0);
    VecNui m_NCells         = VecNui(0);
    VecNui m_NNodes         = VecNui(0);
    UInt   m_NTotalCells    = 1u;
    UInt   m_NTotalNodes    = 1u;
    Real_t m_CellSize       = Real_t(1);
    Real_t m_InvCellSize    = Real_t(1);
    Real_t m_HalfCellSize   = Real_t(0.5);
    Real_t m_CellSizeSqr    = Real_t(1);
    Real_t m_InvCellSizeSqr = Real_t(1);
    Real_t m_CellVolume     = Real_t(1);

    StdVT_UInt m_ParticleIdxSortedByCell;
    bool       m_bCellIdxNeedResize = false; // to track and resize the m_CellParticleIdx array

    Array<N, StdVT_UInt>                m_ParticleIdxInCell;
    Array<N, ParallelObjects::SpinLock> m_Lock;
};
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace NTCodeBase
