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
#include <LibCommon/Grid/Grid.h>
#include <LibCommon/Array/Array.h>
#include <LibCommon/Animation/Animation.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace GeometryObjects
{
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Base class
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
class GeometryObject
{
protected:
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
public:
    GeometryObject() = default;
    static constexpr UInt objDimension() noexcept { return static_cast<UInt>(N); }

    virtual String   name() = 0;
    virtual RealType signedDistance(const VecN& ppos0, bool bNegativeInside = true) const = 0;
    virtual VecN     getAABBMin() const;
    virtual VecN     getAABBMax() const;

    VecN gradSignedDistance(const VecN& ppos, bool bNegativeInside = true, RealType dx = RealType(1e-4)) const;
    bool           isInside(const VecN& ppos, bool bNegativeInside = true) const { return signedDistance(ppos, bNegativeInside) < 0; }

    void setTranslation(const VecN& translation);
    void setRotation(const VecX<N + 1, RealType>& rotation);
    void setUniformScale(const RealType scaleVal);
    void resetTransformation();

    auto& getAnimation() { return m_Animations; }
    auto transformed() const { return m_bTransformed; }
    const auto& getTransformationMatrix() const { return m_TransformationMatrix; }

    VecN transform(const VecN& ppos) const;
    VecN invTransform(const VecN& ppos) const;

    virtual bool updateTransformation(UInt frame = 0, RealType frameFraction = RealType(0), RealType frameDuration = RealType(1.0 / 30.0));

    VecN getLinearVelocity(const VecN& ppos) const;

protected:
    virtual void parseParameters(const JParams& jParams);
    void         updateIntrinsicTransformation();
    ////////////////////////////////////////////////////////////////////////////////
    bool m_bTransformed = false;

    RealType m_LastTime    = 0;
    RealType m_CurrentTime = 0;

    // intrinsic object transformation
    VecN       m_IntrinsicTranslation          = VecN(0);
    VecNp1     m_IntrinsicRotation             = VecNp1(VecN(1), 0);
    RealType   m_UniformScale                  = RealType(1.0);
    MatNp1xNp1 m_IntrinsicTransformationMatrix = MatNp1xNp1(1.0);

    // current transformation = animationTransformation(t) * intrinsicTransformation
    MatNp1xNp1 m_TransformationMatrix    = MatNp1xNp1(1.0);
    MatNp1xNp1 m_InvTransformationMatrix = MatNp1xNp1(1.0);

    StdVT<RigidBodyAnimation<N, RealType>> m_Animations;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Objects
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
class BoxObject : public GeometryObject<N, RealType>
{
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
public:
    BoxObject() = delete;
    BoxObject(const JParams& jParams) { parseParameters(jParams); }

    virtual String   name() override { return String("BoxObject"); }
    virtual RealType signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;

    virtual VecN getAABBMin() const override;
    virtual VecN getAABBMax() const override;

    void        setOriginalBox(const VecN& bMin, const VecN& bMax) { m_BoxMin = bMin; m_BoxMax = bMax; }
    const auto& originalBoxMin() const { return m_BoxMin; }
    const auto& originalBoxMax() const { return m_BoxMax; }

protected:
    virtual void parseParameters(const JParams& jParams) override;

    VecN m_BoxMin = VecN(-1.0);
    VecN m_BoxMax = VecN(1.0);
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
class SphereObject : public GeometryObject<N, RealType>
{
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
public:
    SphereObject() = delete;
    SphereObject(const JParams& jParams) { parseParameters(jParams); }
    virtual String   name() override { if(N == 2) { return String("CircleObject"); } else { return String("SphereObject"); } }
    virtual RealType signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
class TorusObject : public GeometryObject<N, RealType>
{
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
public:
    TorusObject() = delete;
    TorusObject(const JParams& jParams) { parseParameters(jParams); }
    virtual String   name() override { return String("TorusObject"); }
    virtual RealType signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
    void setRingRadius(RealType ringRadius) { __NT_REQUIRE(ringRadius > 0); m_RingRadius = ringRadius; }
protected:
    virtual void parseParameters(const JParams& jParams) override;
    RealType m_OuterRadius = RealType(1);
    RealType m_RingRadius  = RealType(0.25);
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
class Torus28Object : public TorusObject<N, RealType>
{
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
public:
    Torus28Object() = delete;
    Torus28Object(const JParams& jParams) : TorusObject<N, RealType>(jParams) { parseParameters(jParams); }
    virtual RealType signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
class Torus2InfObject : public TorusObject<N, RealType>
{
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
public:
    Torus2InfObject() = delete;
    Torus2InfObject(const JParams& jParams) : TorusObject<N, RealType>(jParams) { parseParameters(jParams); }
    virtual RealType signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
class Torus88Object : public TorusObject<N, RealType>
{
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
public:
    Torus88Object() = delete;
    Torus88Object(const JParams& jParams) : TorusObject<N, RealType>(jParams) { parseParameters(jParams); }
    virtual RealType signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
class TorusInfInfObject : public TorusObject<N, RealType>
{
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
public:
    TorusInfInfObject() = delete;
    TorusInfInfObject(const JParams& jParams) : TorusObject<N, RealType>(jParams) { parseParameters(jParams); }
    virtual RealType signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
class CylinderObject : public GeometryObject<N, RealType>
{
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
public:
    CylinderObject() = delete;
    CylinderObject(const JParams& jParams) { parseParameters(jParams); }
    virtual String   name() override { return String("CylinderObject"); }
    virtual RealType signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
    void setRadius(RealType radius) { m_Radius = radius; }
protected:
    virtual void parseParameters(const JParams& jParams) override;
    RealType m_Radius = RealType(0.5);
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
class ConeObject : public GeometryObject<N, RealType>
{
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
public:
    ConeObject() = delete;
    ConeObject(const JParams& jParams) { parseParameters(jParams); }
    virtual String   name() override { return String("ConeObject"); }
    virtual RealType signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
    void setBaseRadius(RealType radius) { m_Radius = radius; }
protected:
    virtual void parseParameters(const JParams& jParams) override;
    RealType m_Radius = RealType(1);
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
class PlaneObject : public GeometryObject<N, RealType>
{
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
public:
    PlaneObject() = delete;
    PlaneObject(const JParams& jParams) { parseParameters(jParams); }
    virtual String   name() override { return String("PlaneObject"); }
    virtual RealType signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;

    void setNormal(const VecN& normal) { m_Normal = normal; }
    void setOffset(RealType offset) { m_Offset = offset; }
protected:
    virtual void parseParameters(const JParams& jParams) override;
    VecN     m_Normal = VecN(0);
    RealType m_Offset = RealType(0);
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
class TriangleObject : public GeometryObject<N, RealType>
{
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
public:
    TriangleObject() = delete;
    TriangleObject(const JParams& jParams) { parseParameters(jParams); }
    virtual String   name() override { return String("TriangleObject"); }
    virtual RealType signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
    template<class IndexType> void setVertex(IndexType idx, const VecN& vertex) { m_Vertices[idx] = vertex; }
protected:
    virtual void parseParameters(const JParams& jParams) override;
    VecN m_Vertices[3];
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
class HexagonObject : public GeometryObject<N, RealType>
{
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
public:
    HexagonObject() = delete;
    HexagonObject(const JParams& jParams) { parseParameters(jParams); }
    virtual String   name() override { return String("HexagonObject"); }
    virtual RealType signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
class TriangularPrismObject : public GeometryObject<N, RealType>
{
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
public:
    TriangularPrismObject() = delete;
    TriangularPrismObject(const JParams& jParams) { parseParameters(jParams); }
    virtual String   name() override { return String("TriangularPrismObject"); }
    virtual RealType signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
    void setWidth(RealType ratio) { m_Width = ratio; }
protected:
    virtual void parseParameters(const JParams& jParams) override;
    RealType m_Width = RealType(0.5);
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
class HexagonalPrismObject : public GeometryObject<N, RealType>
{
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
public:
    HexagonalPrismObject() = delete;
    HexagonalPrismObject(const JParams& jParams) { parseParameters(jParams); }
    virtual String   name() override { return String("HexagonalPrismObject"); }
    virtual RealType signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
    void             setWidth(RealType ratio) { m_Width = ratio; }
protected:
    virtual void parseParameters(const JParams& jParams) override;
    RealType m_Width = RealType(0.5);
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
class CapsuleObject : public GeometryObject<N, RealType>
{
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
public:
    CapsuleObject() = delete;
    CapsuleObject(const JParams& jParams) { parseParameters(jParams); }
    virtual String   name() override { return String("CapsuleObject"); }
    virtual RealType signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
    void setRadius(RealType radius) { m_Radius = radius; }
protected:
    virtual void parseParameters(const JParams& jParams) override;
    RealType m_Radius = RealType(0.5);
    VecN     m_Start  = VecN(0);
    VecN     m_End    = VecN(0);
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
class EllipsoidObject : public GeometryObject<N, RealType>
{
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
public:
    EllipsoidObject() = delete;
    EllipsoidObject(const JParams& jParams) { parseParameters(jParams); }
    virtual String   name() override { return String("SphereObject"); }
    virtual RealType signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
    void setRadiusRatio(const VecN& scale) { m_RadiusRatio = scale; }
protected:
    virtual void parseParameters(const JParams& jParams) override;
    VecN m_RadiusRatio = VecN(1.0);
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
class TriMeshObject : public GeometryObject<N, RealType>
{
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
public:
    TriMeshObject() = delete;
    TriMeshObject(const JParams& jParams) { parseParameters(jParams); }
    virtual String   name() override { return String("TriMeshObject"); }
    virtual RealType signedDistance(const VecN&, bool bNegativeInside = true) const override;

    auto& meshFile() { return m_TriMeshFile; }
    auto& sdfStep() { return m_Step; }
    void computeSDF();
protected:
    virtual void parseParameters(const JParams& jParams) override;
    bool               m_bSDFGenerated = false;
    String             m_TriMeshFile   = String("");
    RealType           m_Step          = RealType(1.0 / 256.0);
    Array<N, RealType> m_SDFData;
    Grid<N, RealType>  m_Grid3D;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
enum CSGOperations
{
    Overwrite,
    Union,
    Subtraction,
    Intersection,
    BlendExp,
    BlendPoly
};

enum DomainDeformation
{
    None,
    Twist,
    CheapBend
};

template<Int N, class RealType>
class CSGObject : public GeometryObject<N, RealType>
{
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
public:
    struct CSGData
    {
        SharedPtr<GeometryObject<N, RealType>> obj = nullptr;
        CSGOperations                          op  = Union;
    };
    CSGObject() = delete;
    CSGObject(const JParams& jParams) { parseParameters(jParams); }
    virtual String   name() override { return String("CSGObject"); }
    virtual RealType signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;

    void addObject(const CSGData& obj) { m_Objects.push_back(obj); }
    void addObject(const SharedPtr<GeometryObject<N, RealType>>& obj, CSGOperations op = Union) { addObject({ obj, op }); }
    void setDeformOp(DomainDeformation deformOp) { m_DeformOp = deformOp; }

protected:
    VecN domainDeform(const VecN& ppos) const;
    VecN twist(const VecN& ppos) const;
    VecN cheapBend(const VecN& ppos) const;

    ////////////////////////////////////////////////////////////////////////////////
    StdVT<CSGData>    m_Objects;
    DomainDeformation m_DeformOp = None;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
}   // end namespace GeometryObjects
