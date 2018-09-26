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

#include <LibCommon/Grid/Grid.h>
#include <LibCommon/ParallelHelpers/Scheduler.h>
#include <LibCommon/Utils/MathHelpers.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Grid<N, RealType>::setGrid(const VecN& bMin, const VecN& bMax, RealType cellSize, RealType clampEdge /*= 0*/)
{
    m_BMin        = ToVecNa(bMin);
    m_BMax        = ToVecNa(bMax);
    m_ClampedBMin = ToVecNa(bMin) + clampEdge + MEpsilon<RealType>();
    m_ClampedBMax = ToVecNa(bMax) - (clampEdge + MEpsilon<RealType>());
    setCellSize(cellSize);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Grid<N, RealType>::setCellSize(RealType cellSize)
{
    assert(cellSize > 0);
    m_CellSize       = cellSize;
    m_InvCellSize    = RealType(1.0) / m_CellSize;
    m_HalfCellSize   = RealType(0.5) * m_CellSize;
    m_CellSizeSqr    = m_CellSize * m_CellSize;
    m_InvCellSizeSqr = RealType(1.0) / m_CellSizeSqr;

    m_CellVolume  = MathHelpers::pow(m_CellSize, N);
    m_NTotalCells = 1;

    for(Int i = 0; i < m_NCells.length(); ++i) {
        m_NCells[i] = static_cast<UInt>(ceil((m_BMax[i] - m_BMin[i]) / m_CellSize));
        m_NNodes[i] = m_NCells[i] + 1u;

        m_NTotalCells *= m_NCells[i];
        m_NTotalNodes *= m_NNodes[i];
    }

    m_bCellIdxNeedResize = true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Grid<N, RealType>::getGridCoordinate(const StdVT_VecNa& positions, StdVT_VecNa& gridCoordinates) const
{
    assert(positions.size() == gridCoordinates.size());
    Scheduler::parallel_for(positions.size(), [&](size_t p) { gridCoordinates[p] = getGridCoordinate(positions[p]); });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
inline bool Grid<N, RealType>::isInsideGrid(const VecNa& ppos) const noexcept
{
    for(Int d = 0; d < N; ++d) {
        if(ppos[d] < m_BMin[d] || ppos[d] > m_BMax[d]) {
            return false;
        }
    }
    return true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
inline bool Grid<N, RealType>::isInsideClampedBoundary(const VecNa& ppos) const noexcept
{
    for(Int d = 0; d < N; ++d) {
        if(ppos[d] < m_ClampedBMin[d] || ppos[d] > m_ClampedBMax[d]) {
            return false;
        }
    }
    return true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
inline VecX<N, RealType> Grid<N, RealType>::constrainedBoundaryPosition(const VecNa& position) const noexcept
{
    auto constrainedPos = position;
    for(Int d = 0; d < N; ++d) {
        if(constrainedPos[d] < m_BMin[d]) {
            constrainedPos[d] = m_BMin[d];
        } else if(constrainedPos[d] > m_BMax[d]) {
            constrainedPos[d] = m_BMax[d];
        }
    }
    return constrainedPos;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
inline VecX<N, RealType> Grid<N, RealType>::constrainedClampedBoundaryPosition(const VecNa& position) const noexcept
{
    auto constrainedPos = position;
    for(Int d = 0; d < N; ++d) {
        if(constrainedPos[d] < m_ClampedBMin[d]) {
            constrainedPos[d] = m_ClampedBMin[d];
        } else if(constrainedPos[d] > m_ClampedBMax[d]) {
            constrainedPos[d] = m_ClampedBMax[d];
        }
    }
    return constrainedPos;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
inline void Grid<N, RealType>::constrainedBoundaryInPlace(VecNa& position) const noexcept
{
    for(Int d = 0; d < N; ++d) {
        if(position[d] < m_BMin[d]) {
            position[d] = m_BMin[d];
        } else if(position[d] > m_BMax[d]) {
            position[d] = m_BMax[d];
        }
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
inline void Grid<N, RealType>::constrainedClampedBoundaryInPlace(VecNa& position) const noexcept
{
    for(Int d = 0; d < N; ++d) {
        if(position[d] < m_ClampedBMin[d]) {
            position[d] = m_ClampedBMin[d];
        } else if(position[d] > m_ClampedBMax[d]) {
            position[d] = m_ClampedBMax[d];
        }
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Grid<N, RealType>::constrainToGridBoundary(StdVT_VecNa& positions)
{
    Scheduler::parallel_for(positions.size(),
                            [&](size_t p) {
                                auto pos   = positions[p];
                                auto dirty = false;
                                for(Int d = 0; d < N; ++d) {
                                    if(pos[d] < m_BMin[d]) {
                                        dirty  = true;
                                        pos[d] = m_BMin[d];
                                    } else if(pos[d] > m_BMax[d]) {
                                        dirty  = true;
                                        pos[d] = m_BMax[d];
                                    }
                                }
                                if(dirty) {
                                    positions[p] = pos;
                                }
                            });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Grid<N, RealType>::constrainToClampedBoundary(StdVT_VecNa& positions)
{
    Scheduler::parallel_for(positions.size(),
                            [&](size_t p) {
                                auto pos   = positions[p];
                                auto dirty = false;
                                for(Int d = 0; d < N; ++d) {
                                    if(pos[d] < m_ClampedBMin[d]) {
                                        dirty  = true;
                                        pos[d] = m_ClampedBMin[d];
                                    } else if(pos[d] > m_ClampedBMax[d]) {
                                        dirty  = true;
                                        pos[d] = m_ClampedBMax[d];
                                    }
                                }
                                if(dirty) {
                                    positions[p] = pos;
                                }
                            });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Grid<N, RealType>::collectIndexToCells(const StdVT_VecNa& positions)
{
    if(m_bCellIdxNeedResize) {
        m_ParticleIdxInCell.resize(getNCells());
        m_Lock.resize(getNCells());
        m_bCellIdxNeedResize = false;
    }

    for(auto& cell : m_ParticleIdxInCell.flatData()) {
        cell.resize(0);
    }

    Scheduler::parallel_for(static_cast<UInt>(positions.size()),
                            [&](UInt p) {
                                auto cellIdx = getCellIdx<Int>(positions[p]);
                                m_Lock(cellIdx).lock();
                                m_ParticleIdxInCell(cellIdx).push_back(p);
                                m_Lock(cellIdx).unlock();
                            });

    ////////////////////////////////////////////////////////////////////////////////
    // reset particle index vector
    m_ParticleIdxSortedByCell.resize(0);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Grid<N, RealType>::collectIndexToCells(const StdVT_VecNa& positions, StdVT<VecX<N, Int>>& particleCellIdx)
{
    assert(positions.size() == particleCellIdx.size());
    if(m_bCellIdxNeedResize) {
        m_ParticleIdxInCell.resize(getNCells());
        m_Lock.resize(getNCells());
        m_bCellIdxNeedResize = false;
    }

    for(auto& cell : m_ParticleIdxInCell.flatData()) {
        cell.resize(0);
    }

    Scheduler::parallel_for(static_cast<UInt>(positions.size()),
                            [&](UInt p) {
                                auto cellIdx       = getCellIdx<Int>(positions[p]);
                                particleCellIdx[p] = cellIdx;

                                m_Lock(cellIdx).lock();
                                m_ParticleIdxInCell(cellIdx).push_back(p);
                                m_Lock(cellIdx).unlock();
                            });

    ////////////////////////////////////////////////////////////////////////////////
    // reset particle index vector
    m_ParticleIdxSortedByCell.resize(0);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Grid<N, RealType>::collectIndexToCells(const StdVT_VecNa& positions, StdVT_VecNa& gridCoordinates)
{
    assert(positions.size() == gridCoordinates.size());
    if(m_bCellIdxNeedResize) {
        m_ParticleIdxInCell.resize(getNCells());
        m_Lock.resize(getNCells());
        m_bCellIdxNeedResize = false;
    }

    for(auto& cell : m_ParticleIdxInCell.flatData()) {
        cell.resize(0);
    }

    Scheduler::parallel_for(static_cast<UInt>(positions.size()),
                            [&](UInt p) {
                                auto cellPos       = getGridCoordinate(positions[p]);
                                auto cellIdx       = VecX<N, Int>(cellPos);
                                gridCoordinates[p] = cellPos;

                                m_Lock(cellIdx).lock();
                                m_ParticleIdxInCell(cellIdx).push_back(p);
                                m_Lock(cellIdx).unlock();
                            });

    ////////////////////////////////////////////////////////////////////////////////
    // reset particle index vector
    m_ParticleIdxSortedByCell.resize(0);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Grid<N, RealType>::getNeighborList(const StdVT_VecNa& positions, StdVT<StdVT_UInt>& neighborList, Int cellSpan /*= 1*/)
{
    Scheduler::parallel_for(positions.size(), [&](size_t p) { getNeighborList(positions[p], neighborList[p], cellSpan); });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Grid<N, RealType>::getNeighborList(const VecNa& ppos, StdVT_UInt& neighborList, Int cellSpan /*= 1*/)
{
    neighborList.resize(0);
    if constexpr(N == 2) {
        Vec2i cellIdx = getCellIdx<Int>(ppos);
        for(Int lj = -cellSpan; lj <= cellSpan; ++lj) {
            for(Int li = -cellSpan; li <= cellSpan; ++li) {
                const auto neighborCellIdx = cellIdx + Vec2i(li, lj);

                if(!isValidCell(neighborCellIdx)) {
                    continue;
                }

                const auto& cell = m_ParticleIdxInCell(neighborCellIdx);
                if(cell.size() > 0) {
                    neighborList.insert(neighborList.end(), cell.begin(), cell.end());
                }
            }
        }               // end loop over neighbor cells
    } else {
        Vec3i cellIdx = getCellIdx<Int>(ppos);
        for(Int lk = -cellSpan; lk <= cellSpan; ++lk) {
            for(Int lj = -cellSpan; lj <= cellSpan; ++lj) {
                for(Int li = -cellSpan; li <= cellSpan; ++li) {
                    const auto neighborCellIdx = cellIdx + Vec3i(li, lj, lk);

                    if(!isValidCell(neighborCellIdx)) {
                        continue;
                    }

                    const auto& cell = m_ParticleIdxInCell(neighborCellIdx);
                    if(cell.size() > 0) {
                        neighborList.insert(neighborList.end(), cell.begin(), cell.end());
                    }
                }
            }
        }               // end loop over neighbor cells
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Grid<N, RealType>::getNeighborList(const StdVT_VecNa& positions, StdVT<StdVT_UInt>& neighborList, RealType d2, Int cellSpan /*= 1*/)
{
    Scheduler::parallel_for(positions.size(), [&](size_t p) { getNeighborList(positions, positions[p], neighborList[p], d2, cellSpan); });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Grid<N, RealType>::getNeighborList(const StdVT_VecNa& positions, const VecNa& ppos, StdVT_UInt& neighborList, RealType d2, Int cellSpan /*= 1*/)
{
    neighborList.resize(0);

    if constexpr(N == 2) {
        Vec2i cellIdx = getCellIdx<Int>(ppos);
        for(Int lj = -cellSpan; lj <= cellSpan; ++lj) {
            for(Int li = -cellSpan; li <= cellSpan; ++li) {
                const auto neighborCellIdx = cellIdx + Vec2i(li, lj);

                if(!isValidCell(neighborCellIdx)) {
                    continue;
                }

                const auto& cell = m_ParticleIdxInCell(neighborCellIdx);
                if(cell.size() > 0) {
                    for(UInt q : cell) {
                        const auto pqd2 = glm::length2(ppos - positions[q]);
                        if(pqd2 > 0 && pqd2 < d2) {
                            neighborList.push_back(q);
                        }
                    }
                }
            }
        }               // end loop over neighbor cells
    } else {
        Vec3i cellIdx = getCellIdx<Int>(ppos);
        for(Int lk = -cellSpan; lk <= cellSpan; ++lk) {
            for(Int lj = -cellSpan; lj <= cellSpan; ++lj) {
                for(Int li = -cellSpan; li <= cellSpan; ++li) {
                    const auto neighborCellIdx = cellIdx + Vec3i(li, lj, lk);

                    if(!isValidCell(neighborCellIdx)) {
                        continue;
                    }

                    const auto& cell = m_ParticleIdxInCell(neighborCellIdx);

                    if(cell.size() > 0) {
                        for(UInt q : cell) {
                            const auto pqd2 = glm::length2(ppos - positions[q]);
                            if(pqd2 > 0 && pqd2 < d2) {
                                neighborList.push_back(q);
                            }
                        }
                    }
                }
            }
        }               // end loop over neighbor cells
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Grid<N, RealType>::sortData(StdVT_VecNa& data)
{
    const auto& sortedIdx = getParticleIdxSortedByCell();
    assert(sortedIdx.size() == data.size());

    StdVT_VecN tmp(data.begin(), data.end());
    std::transform(sortedIdx.cbegin(), sortedIdx.cend(),
#ifdef _MSC_VER
                   stdext::unchecked_array_iterator<VecNa*>(data.data()),
#else
                   data.data(),
#endif
                   [&](UInt i) { return tmp[i]; });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
const StdVT_UInt& Grid<N, RealType>::getParticleIdxSortedByCell()
{
    if(m_ParticleIdxSortedByCell.size() > 0) {
        return m_ParticleIdxSortedByCell;
    }
    for(auto& cell : m_ParticleIdxInCell.flatData()) {
        if(cell.size() > 0) {
            m_ParticleIdxSortedByCell.insert(m_ParticleIdxSortedByCell.end(), cell.begin(), cell.end());
        }
    }
    return m_ParticleIdxSortedByCell;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(Grid)
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
