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
#include <utility>
#include <algorithm>
#include <functional>

#if !defined(Q_MOC_RUN)
#include <tbb/tbb.h>
#endif

//#define NT_NO_PARALLEL

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NTCodeBase::ParallelExec {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
inline void warmUp() {
    tbb::parallel_for(tbb::blocked_range<Int>(0, 1048576),
                      [&](const tbb::blocked_range<Int>& r) {
                          for(Int i = r.begin(), iEnd = r.end(); i != iEnd; ++i) {
                              volatile int x = 0;
                              (void)x;
                          }
                      });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class IndexType, class Function>
void run(IndexType beginIdx, IndexType endIdx, Function&& function) {
#if defined(NT_NO_PARALLEL) || defined(NT_DISABLE_PARALLEL)
    for(IndexType i = beginIdx; i < endIdx; ++i) {
        function(i);
    }
#else
    tbb::parallel_for(tbb::blocked_range<IndexType>(beginIdx, endIdx),
                      [&](const tbb::blocked_range<IndexType>& r) {
                          for(IndexType i = r.begin(), iEnd = r.end(); i < iEnd; ++i) {
                              function(i);
                          }
                      });
#endif
}

template<class IndexType, class Function>
void run(IndexType endIdx, Function&& function) {
    ParallelExec::run(IndexType(0), endIdx, std::forward<Function>(function));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// parallel for 2D
template<class IndexType, class Function>
void run_row_major(IndexType beginIdxX, IndexType endIdxX,
                   IndexType beginIdxY, IndexType endIdxY,
                   Function&& function) {
    ParallelExec::run(beginIdxX, endIdxX,
                      [&](IndexType i) {
                          for(IndexType j = beginIdxY; j < endIdxY; ++j) {
                              function(i, j);
                          }
                      });
}

template<class IndexType, class Function>
void run(IndexType beginIdxX, IndexType endIdxX,
         IndexType beginIdxY, IndexType endIdxY,
         Function&& function) {
    ParallelExec::run(beginIdxY, endIdxY,
                      [&](IndexType j) {
                          for(IndexType i = beginIdxX; i < endIdxX; ++i) {
                              function(i, j);
                          }
                      });
}

template<class IndexType, class Function>
void run_row_major(const Vec2<IndexType>& endIdx, Function&& function) {
    ParallelExec::run_row_major(IndexType(0), endIdx[0], IndexType(0), endIdx[1], std::forward<Function>(function));
}

template<class IndexType, class Function>
void run(const Vec2<IndexType>& endIdx, Function&& function) {
    ParallelExec::run(IndexType(0), endIdx[0], IndexType(0), endIdx[1], std::forward<Function>(function));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// parallel for 3D
template<class IndexType, class Function>
void run_row_major(IndexType beginIdxX, IndexType endIdxX,
                   IndexType beginIdxY, IndexType endIdxY,
                   IndexType beginIdxZ, IndexType endIdxZ,
                   Function&& function) {
    ParallelExec::run(beginIdxX, endIdxX,
                      [&](IndexType i) {
                          for(IndexType j = beginIdxY; j < endIdxY; ++j) {
                              for(IndexType k = beginIdxZ; k < endIdxZ; ++k) {
                                  function(i, j, k);
                              }
                          }
                      });
}

template<class IndexType, class Function>
void run(IndexType beginIdxX, IndexType endIdxX,
         IndexType beginIdxY, IndexType endIdxY,
         IndexType beginIdxZ, IndexType endIdxZ,
         Function&& function) {
    ParallelExec::run(beginIdxZ, endIdxZ,
                      [&](IndexType k) {
                          for(IndexType j = beginIdxY; j < endIdxY; ++j) {
                              for(IndexType i = beginIdxX; i < endIdxX; ++i) {
                                  function(i, j, k);
                              }
                          }
                      });
}

template<class IndexType, class Function>
void run_row_major(const Vec3<IndexType>& endIdx, Function&& function) {
    ParallelExec::run_row_major(IndexType(0), endIdx[0], IndexType(0), endIdx[1], IndexType(0), endIdx[2], std::forward<Function>(function));
}

template<class IndexType, class Function>
void run(const Vec3<IndexType>& endIdx, Function&& function) {
    ParallelExec::run(IndexType(0), endIdx[0], IndexType(0), endIdx[1], IndexType(0), endIdx[2], std::forward<Function>(function));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // namespace NTCodeBase::ParallelExec
