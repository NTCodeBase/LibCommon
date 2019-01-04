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

#include <tbb/tbb.h>

#include <utility>
#include <algorithm>
#include <functional>

//#define __NT_NO_PARALLEL
//#define __NT_PARALLEL_FOR_UNROLL

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NTCodeBase::Scheduler {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
inline void warmUp() {
    tbb::parallel_for(tbb::blocked_range<Int>(0, 1048576),
                      [&](const tbb::blocked_range<Int>& r) {
                          for(Int i = r.begin(), iEnd = r.end(); i != iEnd; ++i) {
                              volatile int x;
                              (void)x;
                          }
                      });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class IndexType, class Function>
inline void parallel_for(IndexType beginIdx, IndexType endIdx, Function&& function) {
#if defined(__NT_NO_PARALLEL) || defined(__NT_DISABLE_PARALLEL)
    for(IndexType i = beginIdx; i < endIdx; ++i) {
        function(i);
    }
#else
#  ifdef __NT_PARALLEL_FOR_UNROLL
    tbb::parallel_for(tbb::blocked_range<IndexType>(beginIdx, endIdx),
                      [&](const tbb::blocked_range<IndexType>& r) {
                          auto i = r.begin();
                          auto n = (r.end() - r.begin()) / IndexType(4);
                          for(IndexType j = 0; j < n; i += IndexType(4), ++j) {
                              function(i);
                              function(i + IndexType(1));
                              function(i + IndexType(2));
                              function(i + IndexType(3));
                          }
                          for(auto iEnd = r.end(); i < iEnd; ++i) {
                              function(i);
                          }
                      });
#  else
    tbb::parallel_for(tbb::blocked_range<IndexType>(beginIdx, endIdx),
                      [&](const tbb::blocked_range<IndexType>& r) {
                          for(IndexType i = r.begin(), iEnd = r.end(); i < iEnd; ++i) {
                              function(i);
                          }
                      });
#  endif
#endif
}

template<class IndexType, class Function>
inline void parallel_for(IndexType endIdx, Function&& function) {
    Scheduler::parallel_for(IndexType(0), endIdx, std::forward<Function>(function));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// parallel for 2D
template<class IndexType, class Function>
inline void parallel_for_row_major(IndexType beginIdxX, IndexType endIdxX,
                                   IndexType beginIdxY, IndexType endIdxY,
                                   Function&& function) {
    Scheduler::parallel_for(beginIdxX, endIdxX,
                            [&](IndexType i) {
                                for(IndexType j = beginIdxY; j < endIdxY; ++j) {
                                    function(i, j);
                                }
                            });
}

template<class IndexType, class Function>
inline void parallel_for(IndexType beginIdxX, IndexType endIdxX,
                         IndexType beginIdxY, IndexType endIdxY,
                         Function&& function) {
    Scheduler::parallel_for(beginIdxY, endIdxY,
                            [&](IndexType j) {
                                for(IndexType i = beginIdxX; i < endIdxX; ++i) {
                                    function(i, j);
                                }
                            });
}

template<class IndexType, class Function>
inline void parallel_for_row_major(const Vec2<IndexType>& endIdx, Function&& function) {
    Scheduler::parallel_for_row_major(IndexType(0), endIdx[0], IndexType(0), endIdx[1], std::forward<Function>(function));
}

template<class IndexType, class Function>
inline void parallel_for(const Vec2<IndexType>& endIdx, Function&& function) {
    Scheduler::parallel_for(IndexType(0), endIdx[0], IndexType(0), endIdx[1], std::forward<Function>(function));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// parallel for 3D
template<class IndexType, class Function>
inline void parallel_for_row_major(IndexType beginIdxX, IndexType endIdxX,
                                   IndexType beginIdxY, IndexType endIdxY,
                                   IndexType beginIdxZ, IndexType endIdxZ,
                                   Function&& function) {
    Scheduler::parallel_for(beginIdxX, endIdxX,
                            [&](IndexType i) {
                                for(IndexType j = beginIdxY; j < endIdxY; ++j) {
                                    for(IndexType k = beginIdxZ; k < endIdxZ; ++k) {
                                        function(i, j, k);
                                    }
                                }
                            });
}

template<class IndexType, class Function>
inline void parallel_for(IndexType beginIdxX, IndexType endIdxX,
                         IndexType beginIdxY, IndexType endIdxY,
                         IndexType beginIdxZ, IndexType endIdxZ,
                         Function&& function) {
    Scheduler::parallel_for(beginIdxZ, endIdxZ,
                            [&](IndexType k) {
                                for(IndexType j = beginIdxY; j < endIdxY; ++j) {
                                    for(IndexType i = beginIdxX; i < endIdxX; ++i) {
                                        function(i, j, k);
                                    }
                                }
                            });
}

template<class IndexType, class Function>
inline void parallel_for_row_major(const Vec3<IndexType>& endIdx, Function&& function) {
    Scheduler::parallel_for_row_major(IndexType(0), endIdx[0], IndexType(0), endIdx[1], IndexType(0), endIdx[2], std::forward<Function>(function));
}

template<class IndexType, class Function>
inline void parallel_for(const Vec3<IndexType>& endIdx, Function&& function) {
    Scheduler::parallel_for(IndexType(0), endIdx[0], IndexType(0), endIdx[1], IndexType(0), endIdx[2], std::forward<Function>(function));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // namespace NTCodeBase::Scheduler
