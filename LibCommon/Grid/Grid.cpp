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

#include <Grid/Grid.h>
#include <Utils/MathHelpers.h>
#include <ParallelHelpers/Scheduler.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Grid<N, RealType>::setGrid(const VecN& bMin, const VecN& bMax, RealType cellSize)
{
    m_BMin = bMin;
    m_BMax = bMax;
    setCellSize(cellSize);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Grid<N, RealType>::setCellSize(RealType cellSize)
{
    assert(cellSize > 0);
    m_CellSize     = cellSize;
    m_InvCellSize  = RealType(1.0) / m_CellSize;
    m_HalfCellSize = RealType(0.5) * m_CellSize;
    m_CellSizeSqr  = m_CellSize * m_CellSize;

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
void Grid<N, RealType>::getGridCoordinate(const StdVT_VecN& positions, StdVT_VecN& gridCoordinates) const
{
    assert(positions.size() == gridCoordinates.size());
    Scheduler::parallel_for(positions.size(), [&](size_t p) { gridCoordinates[p] = getGridCoordinate(positions[p]); });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
bool Grid<N, RealType>::isInsideGrid(const VecN& ppos) const noexcept
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
void Grid<N, RealType>::constrainToGrid(StdVT_VecN& positions)
{
    Scheduler::parallel_for(positions.size(),
                            [&](size_t p)
                            {
                                VecN pos   = positions[p];
                                bool dirty = false;

                                for(Int d = 0; d < N; ++d) {
                                    if(pos[d] < m_BMin[d] || pos[d] > m_BMax[d]) {
                                        dirty  = true;
                                        pos[d] = MathHelpers::clamp(pos[d], m_BMin[d], m_BMax[d]);
                                    }
                                }

                                if(dirty) {
                                    positions[p] = pos;
                                }
                            });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Grid<N, RealType>::collectIndexToCells(const StdVT_VecN& positions)
{
    if(m_bCellIdxNeedResize) {
        m_ParticleIdxInCell.resize(getNCells());
        m_Lock.resize(getNCells());
        m_bCellIdxNeedResize = false;
    }

    for(auto& cell : m_ParticleIdxInCell.data()) {
        cell.resize(0);
    }

    Scheduler::parallel_for(static_cast<UInt>(positions.size()),
                            [&](UInt p)
                            {
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
void Grid<N, RealType>::collectIndexToCells(const StdVT_VecN& positions, StdVT<VecX<N, Int>>& particleCellIdx)
{
    assert(positions.size() == particleCellIdx.size());
    if(m_bCellIdxNeedResize) {
        m_ParticleIdxInCell.resize(getNCells());
        m_Lock.resize(getNCells());
        m_bCellIdxNeedResize = false;
    }

    for(auto& cell : m_ParticleIdxInCell.data()) {
        cell.resize(0);
    }

    Scheduler::parallel_for(static_cast<UInt>(positions.size()),
                            [&](UInt p)
                            {
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
void Grid<N, RealType>::collectIndexToCells(const StdVT_VecN& positions, StdVT_VecN& gridCoordinates)
{
    assert(positions.size() == gridCoordinates.size());
    if(m_bCellIdxNeedResize) {
        m_ParticleIdxInCell.resize(getNCells());
        m_Lock.resize(getNCells());
        m_bCellIdxNeedResize = false;
    }

    for(auto& cell : m_ParticleIdxInCell.data()) {
        cell.resize(0);
    }

    Scheduler::parallel_for(static_cast<UInt>(positions.size()),
                            [&](UInt p)
                            {
                                auto cellPos       = getCellIdx<RealType>(positions[p]);
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
void Grid<N, RealType>::getNeighborList(const StdVT_VecN& positions, StdVT<StdVT_UInt>& neighborList, Int cellSpan /*= 1*/)
{
    Scheduler::parallel_for(positions.size(), [&](size_t p) { getNeighborList(positions[p], neighborList[p], cellSpan); });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Grid<N, RealType>::getNeighborList(const VecN& ppos, StdVT_UInt& neighborList, Int cellSpan /*= 1*/)
{
    neighborList.resize(0);
    if constexpr(N == 2)
    {
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
void Grid<N, RealType>::getNeighborList(const StdVT_VecN& positions, StdVT<StdVT_UInt>& neighborList, RealType d2, Int cellSpan /*= 1*/)
{
    Scheduler::parallel_for(positions.size(), [&](size_t p) { getNeighborList(positions, positions[p], neighborList[p], d2, cellSpan); });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Grid<N, RealType>::getNeighborList(const StdVT_VecN& positions, const VecN& ppos, StdVT_UInt& neighborList, RealType d2, Int cellSpan /*= 1*/)
{
    neighborList.resize(0);

    if constexpr(N == 2)
    {
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
void Grid<N, RealType>::sortData(StdVT_VecN& data)
{
    const auto& sortedIdx = getParticleIdxSortedByCell();
    assert(sortedIdx.size() == data.size());

    StdVT_VecN tmp(data.begin(), data.end());
    std::transform(sortedIdx.cbegin(), sortedIdx.cend(),
#ifdef _MSC_VER
                   stdext::unchecked_array_iterator<VecN*>(data.data()),
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
    for(auto& cell : m_ParticleIdxInCell.data()) {
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
