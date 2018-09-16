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
#include <LibCommon/ParallelHelpers/ParallelObjects.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// The pool of objects that are used mutually exclusively
// In addition, if no object is available, new object will be created to serve the request
template<class Object>
class DynamicMutexPool
{
public:
    DynamicMutexPool(UInt initSize) : m_Objects(initSize, Object()), m_bAvailable(initSize, 1) {}

    /**
     * \brief Get the next available object idx. If no object is available, insert a new object and return it.
     */
    size_t getAvailableIdx()
    {
        size_t idx = Huge<size_t>();
        m_Lock.lock();
        for(size_t i = 0; i < m_bAvailable.size(); ++i) {
            if(m_bAvailable[i]) {
                m_bAvailable[i] = 0;
                idx = i;
                break;
            }
        }

        // if no object is availble: insert a new object and return it
        if(idx == Huge<size_t>()) {
            idx = m_Objects.size();
            m_bAvailable.push_back(0);
            m_Objects.emplace_back(Object());
        }
        m_Lock.unlock();
        return idx;
    }

    /**
     * \brief Get the object. The parameter idx should be the index of an available object
     */
    template<class IndexType>
    auto& getObj(IndexType idx)
    {
        return m_Objects[idx];
    }

    /**
     * \brief Set object to be available for next use
     */
    template<class IndexType>
    void setAvailable(IndexType idx)
    {
        m_bAvailable[idx] = 1;
    }

private:
    ParallelObjects::SpinLock m_Lock;
    StdVT<Object>             m_Objects;
    StdVT<Int8>               m_bAvailable;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Fixed size pool where each object may be used by more than one threads
template<size_t N, class Object>
class CyclePool
{
public:
    CyclePool() = default;

    /**
     * \brief Get the next object
     */
    Object& getNext()
    {
        m_Lock.lock();
        if(++m_CurrentIdx == N) {
            m_CurrentIdx = 0;
        }
        size_t idx = m_CurrentIdx;
        m_Lock.unlock();
        return m_Objects[idx];
    }

private:
    ParallelObjects::SpinLock m_Lock;
    Object                    m_Objects[N];
    size_t                    m_CurrentIdx = 0;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
