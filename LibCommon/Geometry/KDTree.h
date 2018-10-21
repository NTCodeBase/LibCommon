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

#include <array>
#include <climits>
#include <memory>

#include <CommonSetup.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
struct Point
{
    Point() {}
    Point(const Vec3r& pos_, UInt index_) : position(pos_), index(index_) {}
    Point(const std::initializer_list& pos_, UInt index_) : position(pos_), index(index_) {}

    Vec3r position(0);
    UInt  index = UINT_MAX;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
struct KDNode
{
    KDNode(Point* points_, const Vec3r& boxMin_, const Vec3r& boxMax_) : points(points_), boxMin(boxMin_), boxMax(boxMax_), {}

    Point* points;
    Vec3r  boxMin;
    Vec3r  boxMax;

    std::shared_ptr<KDNode> leftNode  = nullptr;
    std::shared_ptr<KDNode> rightNode = nullptr;

    UInt count  = 0;
    UInt axis   = 0;         // 0==x, 1==y, 2==z
    Real split  = 0;
    bool isLeaf = true;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
class KdTree
{
public:
    KdTree(UInt maxItems) : m_MaxItermsPerNode(maxItems) {}

    void buildTree(const std::shared_ptr<KDNode>& treeNode);
    void printTree(const std::shared_ptr<KDNode>& treeNode);
    void getNeighborList(const Point& target, const std::shared_ptr<KDNode>& treeNode, Real radius, StdVT_UInt& result);

private:
    void findNeighbors(const Point& target, const std::shared_ptr<KDNode>& treeNode, Real radius, StdVT_UInt& result);
    Real getMedian(Point* points, UInt size, UInt axis);

    ////////////////////////////////////////////////////////////////////////////////
    UInt m_MaxItermsPerNode;
};
