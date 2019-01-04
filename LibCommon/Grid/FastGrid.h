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

#include <cassert>
#include <LibCommon/CommonSetup.h>
#include <LibCommon/Math/MathHelpers.h>
#include <LibCommon/Math/FastVec3.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NTCodeBase {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
class FastGrid3 {
public:
    inline FastGrid3() = default;
    ////////////////////////////////////////////////////////////////////////////////
    inline void setGrid(const Vec3<Real_t>& bMin, const Vec3<Real_t>& bMax, Real_t cellSize, Real_t clampEdge = 0) {
        m_BMin        = bMin;
        m_BMax        = bMax;
        m_ClampedBMin = bMin + clampEdge + MEpsilon<Real_t>();
        m_ClampedBMax = bMax - (clampEdge + MEpsilon<Real_t>());
        ////////////////////////////////////////////////////////////////////////////////
        m_CellSize     = cellSize;
        m_InvCellSize  = Real_t(1.0) / m_CellSize;
        m_HalfCellSize = Real_t(0.5) * m_CellSize;
        m_CellSizeSqr  = m_CellSize * m_CellSize;

        m_CellVolume  = MathHelpers::pow(m_CellSize, 3);
        m_NTotalCells = 1;

        for(Int i = 0; i < m_NCells.length(); ++i) {
            m_NCells[i] = static_cast<UInt>(ceil((m_BMax[i] - m_BMin[i]) / m_CellSize));
            m_NNodes[i] = m_NCells[i] + 1u;

            m_NTotalCells *= m_NCells[i];
            m_NTotalNodes *= m_NNodes[i];
        }
    }

    ////////////////////////////////////////////////////////////////////////////////
    inline const auto& getNCells() const noexcept { return m_NCells; }
    inline const auto& getNNodes() const noexcept { return m_NNodes; }
    inline auto        getNTotalCells() const noexcept { return m_NTotalCells; }
    inline auto        getNTotalNodes() const noexcept { return m_NTotalNodes; }
    ////////////////////////////////////////////////////////////////////////////////
    inline const auto& getBMin() const { return m_BMin; }
    inline const auto& getBMax() const { return m_BMax; }
    inline auto getCellSize() const noexcept { return m_CellSize; }
    inline auto getInvCellSize() const noexcept { return m_InvCellSize; }
    inline auto getHalfCellSize() const noexcept { return m_HalfCellSize; }
    inline auto getCellSizeSquared() const noexcept { return m_CellSizeSqr; }
    inline auto getCellVolume() const noexcept { return m_CellVolume; }
    ////////////////////////////////////////////////////////////////////////////////
    inline auto getWorldCoordinate(const Vec3i& nodeIdx) const { return FastVec3<Real_t>(nodeIdx) * m_CellSize + m_BMin; }
    inline auto getWorldCoordinate(int i, int j, int k) const { return FastVec3<Real_t>(i, j, k) * m_CellSize + m_BMin; }
    inline auto getGridCoordinate(const FastVec3<Real_t>& ppos) const { return (ppos - m_BMin) * m_InvCellSize; }
protected:
    FastVec3<Real_t> m_BMin         = FastVec3<Real_t>(-1.0);
    FastVec3<Real_t> m_BMax         = FastVec3<Real_t>(1.0);
    FastVec3<Real_t> m_ClampedBMin  = FastVec3<Real_t>(-1.0);
    FastVec3<Real_t> m_ClampedBMax  = FastVec3<Real_t>(1.0);
    Vec3ui           m_NCells       = Vec3ui(0);
    Vec3ui           m_NNodes       = Vec3ui(0);
    UInt             m_NTotalCells  = 1u;
    UInt             m_NTotalNodes  = 1u;
    Real_t           m_CellSize     = Real_t(1);
    Real_t           m_InvCellSize  = Real_t(1);
    Real_t           m_HalfCellSize = Real_t(0.5);
    Real_t           m_CellSizeSqr  = Real_t(1);
    Real_t           m_CellVolume   = Real_t(1);
};
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace NTCodeBase
