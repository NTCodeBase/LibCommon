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
#include <LibCommon/Array/Array.h>
#include <LibCommon/Utils/MathHelpers.h>
#include <LibCommon/ParallelHelpers/ParallelObjects.h>
#include <cassert>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
class Grid
{
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
public:
    Grid() = default;
    Grid(const VecN& bMin, const VecN& bMax, RealType cellSize) : m_BMin(ToVecNa(bMin)), m_BMax(ToVecNa(bMax)) { setCellSize(cellSize); }

    ////////////////////////////////////////////////////////////////////////////////
    // Setters
    void setGrid(const VecN& bMin, const VecN& bMax, RealType cellSize, RealType clampEdge = 0);
    void setCellSize(RealType cellSize);

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
    IndexType getCellLinearizedIndex(IndexType i, IndexType j) const noexcept
    {
        static_assert(N == 2, "Array dimension != 2");
        return j * static_cast<IndexType>(getNCells()[0]) + i;
    }

    template<class IndexType>
    IndexType getNodeLinearizedIndex(IndexType i, IndexType j) const noexcept
    {
        static_assert(N == 2, "Array dimension != 2");
        return j * static_cast<IndexType>(getNNodes()[0]()) + i;
    }

    template<class IndexType>
    bool isValidCell(IndexType i, IndexType j)  const noexcept
    {
        static_assert(N == 2, "Array dimension != 2");
        return (i >= 0 &&
                j >= 0 &&
                static_cast<UInt>(i) < m_NCells[0] &&
                static_cast<UInt>(j) < m_NCells[1]);
    }

    template<class IndexType>
    bool isValidCell(const Vec2<IndexType>& index) const noexcept
    {
        return isValidCell(index[0], index[1]);
    }

    template<class IndexType>
    bool isValidNode(IndexType i, IndexType j)  const noexcept
    {
        static_assert(N == 2, "Array dimension != 2");
        return (i >= 0 &&
                j >= 0 &&
                static_cast<UInt>(i) < m_NNodes[0] &&
                static_cast<UInt>(j) < m_NNodes[1]);
    }

    template<class IndexType>
    bool isValidNode(const Vec2<IndexType>& index) const noexcept
    {
        return isValidNode(index[0], index[1]);
    }

    // => Grid 2D
    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////
    // Grid 3D =>
    template<class IndexType>
    IndexType getCellLinearizedIndex(IndexType i, IndexType j, IndexType k) const
    {
        static_assert(N == 3, "Array dimension != 3");
        return (k * static_cast<IndexType>(getNCells()[1]) + j) * static_cast<IndexType>(getNCells()[0]) + i;
    }

    template<class IndexType>
    IndexType getNodeLinearizedIndex(IndexType i, IndexType j, IndexType k) const
    {
        static_assert(N == 3, "Array dimension != 3");
        return (k * static_cast<IndexType>(getNNodes()[1]()) + j) * static_cast<IndexType>(getNNodes()[0]) + i;
    }

    template<class IndexType>
    bool isValidCell(IndexType i, IndexType j, IndexType k)  const noexcept
    {
        static_assert(N == 3, "Array dimension != 3");
        return (i >= 0 &&
                j >= 0 &&
                k >= 0 &&
                static_cast<UInt>(i) < m_NCells[0] &&
                static_cast<UInt>(j) < m_NCells[1] &&
                static_cast<UInt>(k) < m_NCells[2]);
    }

    template<class IndexType>
    bool isValidCell(const Vec3<IndexType>& index) const noexcept
    {
        return isValidCell(index[0], index[1], index[2]);
    }

    template<class IndexType>
    bool isValidNode(IndexType i, IndexType j, IndexType k)  const noexcept
    {
        static_assert(N == 3, "Array dimension != 3");
        return (i >= 0 &&
                j >= 0 &&
                k >= 0 &&
                static_cast<UInt>(i) < m_NNodes[0] &&
                static_cast<UInt>(j) < m_NNodes[1] &&
                static_cast<UInt>(k) < m_NNodes[2]);
    }

    template<class IndexType>
    bool isValidNode(const Vec3<IndexType>& index) const noexcept
    {
        return isValidNode(index[0], index[1], index[2]);
    }

    // => Grid 3D
    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////
    template<class IndexType>
    auto getCellIdx(const VecNa& ppos) const noexcept
    {
        VecX<N, IndexType> cellIdx;
        for(Int d = 0; d < N; ++d) {
            cellIdx[d] = static_cast<IndexType>((ppos[d] - m_BMin[d]) / m_CellSize);
        }
        return cellIdx;
    }

    template<class IndexType>
    auto getNearestValidCellIdx(const VecX<N, IndexType>& cellIdx) const noexcept
    {
        VecX<N, IndexType> validCellIdx;
        for(Int d = 0; d < N; ++d) {
            validCellIdx[d] = MathHelpers::clamp<IndexType>(cellIdx[d], 0, static_cast<IndexType>(m_NCells[d]) - 1);
        }
        return validCellIdx;
    }

    template<class IndexType>
    auto getValidCellIdx(const VecNa& ppos) const noexcept
    {
        return getNearestValidCellIdx<IndexType>(getCellIdx<IndexType>(ppos));
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Particle processing
    template<class IndexType>
    auto getWorldCoordinate(const VecX<N, IndexType>& nodeIdx) const
    {
        return VecNa(nodeIdx) * m_CellSize + m_BMin;
    }

    template<class IndexType>
    auto getWorldCoordinate(IndexType i, IndexType j) const
    {
        static_assert(N == 2, "Array dimension != 2");
        return Vec2<RealType>(i, j) * m_CellSize + m_BMin;
    }

    template<class IndexType>
    auto getWorldCoordinate(IndexType i, IndexType j, IndexType k) const
    {
        static_assert(N == 3, "Array dimension != 3");
        return Vec4<RealType>(i, j, k, 0) * m_CellSize + m_BMin;
    }

    auto getGridCoordinate(const VecNa& ppos) const { return (ppos - m_BMin) * m_InvCellSize; }
    void getGridCoordinate(const StdVT_VecNa& positions, StdVT_VecNa& gridCoordinates) const;
    ////////////////////////////////////////////////////////////////////////////////
    inline bool isInsideGrid(const VecNa& ppos) const noexcept;
    inline bool isInsideClampedBoundary(const VecNa& ppos) const noexcept;
    inline VecN constrainedBoundaryPosition(const VecNa& positions) const noexcept;
    inline VecN constrainedClampedBoundaryPosition(const VecNa& positions) const noexcept;
    inline void constrainedBoundaryInPlace(VecNa& positions) const noexcept;
    inline void constrainedClampedBoundaryInPlace(VecNa& positions) const noexcept;

    void constrainToGridBoundary(StdVT_VecNa& positions);
    void constrainToClampedBoundary(StdVT_VecNa& positions);
    void collectIndexToCells(const StdVT_VecNa& positions);
    void collectIndexToCells(const StdVT_VecNa& positions, StdVT<VecX<N, Int>>& particleCellIdx);
    void collectIndexToCells(const StdVT_VecNa& positions, StdVT_VecNa& gridCoordinates);
    void getNeighborList(const StdVT_VecNa& positions, StdVT<StdVT_UInt>& neighborList, Int cellSpan = 1);
    void getNeighborList(const VecNa& ppos, StdVT_UInt& neighborList, Int cellSpan = 1);
    void getNeighborList(const StdVT_VecNa& positions, StdVT<StdVT_UInt>& neighborList, RealType d2, Int cellSpan = 1);
    void getNeighborList(const StdVT_VecNa& positions, const VecNa& ppos, StdVT_UInt& neighborList, RealType d2, Int cellSpan = 1);
    void sortData(StdVT_VecNa& data);

    const StdVT_UInt& getParticleIdxSortedByCell();

    template<class IndexType>
    const auto& getParticleIdxInCell(const VecX<N, IndexType>& cellIdx) const { return m_ParticleIdxInCell(cellIdx); }
protected:
    VecNa    m_BMin           = VecNa(-1.0);
    VecNa    m_BMax           = VecNa(1.0);
    VecNa    m_ClampedBMin    = VecNa(-1.0);
    VecNa    m_ClampedBMax    = VecNa(1.0);
    VecNui   m_NCells         = VecX<N, UInt>(0);
    VecNui   m_NNodes         = VecX<N, UInt>(0);
    UInt     m_NTotalCells    = 1u;
    UInt     m_NTotalNodes    = 1u;
    RealType m_CellSize       = RealType(1);
    RealType m_InvCellSize    = RealType(1);
    RealType m_HalfCellSize   = RealType(0.5);
    RealType m_CellSizeSqr    = RealType(1);
    RealType m_InvCellSizeSqr = RealType(1);
    RealType m_CellVolume     = RealType(1);

    StdVT_UInt m_ParticleIdxSortedByCell;
    bool       m_bCellIdxNeedResize = false; // to track and resize the m_CellParticleIdx array

    Array<N, StdVT_UInt>                m_ParticleIdxInCell;
    Array<N, ParallelObjects::SpinLock> m_Lock;
};
