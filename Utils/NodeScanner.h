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

#include <Setup.h>
#include <Core/Utils/MathHelpers.h>
#include <array>
#include <cassert>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int Span, Int N, class NodeIdxType, class FlatIdxType>
class NodeScanner
{
public:
    NodeScanner(const VecX<N, NodeIdxType>& beginNode = VecX<N, NodeIdxType>(0), FlatIdxType beginFlatIndex = 0) :
        m_BeginNode(beginNode), m_CurrentNode(beginNode), m_CurrentFlatIndex(beginFlatIndex) {}

    /**
     * \brief Get number of scanning nodes
     */
    static constexpr UInt nNodes() { return m_nNodes; }

    /**
     * \brief Check if the scanner reaches the end
     */
    bool notEnd() const { return m_bNotEnd; }
    bool isEnd() const { return !m_bNotEnd; }

    /**
     * \brief Get linearized index
     */
    auto flatIndex() const { return m_CurrentFlatIndex; }

    /**
     * \brief Get current node from [beginNode, beginNode + nNodes())
     */
    const auto& currentNode() const { return m_CurrentNode; }

    /**
     * \brief Get zero-based node
     */
    const auto& basedNode() const { return m_BaseNodes[m_BaseFlatIndex]; }

    /**
     * \brief Increment the linearized index by one
     */
    void operator++()
    {
        ++m_BaseFlatIndex;
        if(m_BaseFlatIndex < m_nNodes) {
            m_CurrentNode = m_BaseNodes[m_BaseFlatIndex] + m_BeginNode;
            ++m_CurrentFlatIndex;
        } else {
            m_bNotEnd = false;
        }
    }

private:
    static constexpr auto populateBaseNodes()
    {
        static_assert(N == 2 || N == 3);
        std::array<VecX<N, NodeIdxType>, MathHelpers::pow(Span, N)> idxList;

        auto idx = UInt(0);
        if constexpr (N == 2) {
            for(NodeIdxType j = 0; j < Span; ++j) {
                for(NodeIdxType i = 0; i < Span; ++i) {
                    idxList[idx++] = VecX<N, NodeIdxType>(i, j);
                }
            }
        } else {
            for(NodeIdxType k = 0; k < Span; ++k) {
                for(NodeIdxType j = 0; j < Span; ++j) {
                    for(NodeIdxType i = 0; i < Span; ++i) {
                        idxList[idx++] = VecX<N, NodeIdxType>(i, j, k);
                    }
                }
            }
        }
        return idxList;
    }

    ////////////////////////////////////////////////////////////////////////////////
    VecX<N, NodeIdxType> m_BeginNode;
    VecX<N, NodeIdxType> m_CurrentNode;
    FlatIdxType          m_CurrentFlatIndex;

    FlatIdxType m_BaseFlatIndex = 0;
    bool        m_bNotEnd       = true;

    // the total number of nodes
    const static inline UInt m_nNodes = static_cast<UInt>(MathHelpers::pow(Span, N));

    // store all based nodes
    const static inline std::array<VecX<N, NodeIdxType>, MathHelpers::pow(Span, N)> m_BaseNodes = populateBaseNodes();
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class NodeIdxType, class FlatIdxType>
auto NodeScanner2(const VecX<N, NodeIdxType>& beginNode = VecX<N, NodeIdxType>(0), FlatIdxType beginFlatIndex = 0)
{
    return NodeScanner<2, N, NodeIdxType, FlatIdxType>(beginNode, beginFlatIndex);
}

template<Int N, class NodeIdxType, class FlatIdxType>
auto NodeScanner3(const VecX<N, NodeIdxType>& beginNode = VecX<N, NodeIdxType>(0), FlatIdxType beginFlatIndex = 0)
{
    return NodeScanner<3, N, NodeIdxType, FlatIdxType>(beginNode, beginFlatIndex);
}

template<Int N, class NodeIdxType, class FlatIdxType>
auto NodeScanner4(const VecX<N, NodeIdxType>& beginNode = VecX<N, NodeIdxType>(0), FlatIdxType beginFlatIndex = 0)
{
    return NodeScanner<4, N, NodeIdxType, FlatIdxType>(beginNode, beginFlatIndex);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
