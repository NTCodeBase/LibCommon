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
#include <LibCommon/Grid/Grid.h>
#include <LibCommon/Array/Array.h>
#include <LibCommon/Animation/Animation.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NTCodeBase {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Base class
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class GeometryObject {
protected:
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIAS
    ////////////////////////////////////////////////////////////////////////////////
public:
    GeometryObject() = default;
    static constexpr UInt objDimension() noexcept { return static_cast<UInt>(N); }

    virtual String name() = 0;
    virtual Real_t signedDistance(const VecN& ppos0, bool bNegativeInside = true) const = 0;
    virtual VecN   getAABBMin() const;
    virtual VecN   getAABBMax() const;

    VecN gradSignedDistance(const VecN& ppos, bool bNegativeInside = true, Real_t dx = Real_t(1e-4)) const;
    bool           isInside(const VecN& ppos, bool bNegativeInside = true) const { return signedDistance(ppos, bNegativeInside) < 0; }

    void setTranslation(const VecN& translation);
    void setRotation(const VecNp1& rotation);
    void setUniformScale(const Real_t scaleVal);
    auto uniformScale() const { return m_UniformScale; }
    void resetTransformation();

    auto& getAnimation() { return m_Animations; }
    auto transformed() const { return m_bTransformed; }
    auto animationTransformed() const { return m_bTransformed && m_Animations.size() > 0; }
    const auto& getTransformationMatrix() const { return m_TransformationMatrix; }
    const auto& getPrevTransformationMatrix() const { return m_PrevTransformationMatrix; }
    const auto& getAnimationTransformationMatrix() const { return m_AnimationTransformationMatrix; }
    const auto& getPrevAnimationTransformationMatrix() const { return m_PrevAnimationTransformationMatrix; }

    bool doneTransformation() const { return m_bDoneTransformation; }
    bool updateTransformation(UInt frame = 0, Real_t frameFraction = Real_t(0));
    VecN transformAnimation(const VecN& ppos) const;
    VecN invTransformAnimation(const VecN& ppos) const;
    VecN transform(const VecN& ppos) const;
    VecN invTransform(const VecN& ppos) const;

protected:
    virtual void parseParameters(const JParams& jParams);
    void         updateIntrinsicTransformation();
    ////////////////////////////////////////////////////////////////////////////////
    bool m_bTransformed        = false;
    bool m_bDoneTransformation = false;

    // intrinsic object transformation
    VecN       m_IntrinsicTranslation          = VecN(0);
    VecNp1     m_IntrinsicRotation             = VecNp1(VecN(1), 0);
    Real_t     m_UniformScale                  = Real_t(1.0);
    MatNp1xNp1 m_IntrinsicTransformationMatrix = MatNp1xNp1(1.0);

    // current transformation = animationTransformation(t) * intrinsicTransformation
    MatNp1xNp1 m_TransformationMatrix     = MatNp1xNp1(1.0);
    MatNp1xNp1 m_InvTransformationMatrix  = MatNp1xNp1(1.0);
    MatNp1xNp1 m_PrevTransformationMatrix = MatNp1xNp1(1.0);

    MatNp1xNp1 m_AnimationTransformationMatrix     = MatNp1xNp1(1.0);
    MatNp1xNp1 m_InvAnimationTransformationMatrix  = MatNp1xNp1(1.0);
    MatNp1xNp1 m_PrevAnimationTransformationMatrix = MatNp1xNp1(1.0);
    StdVT<RigidBodyAnimation<N, Real_t>> m_Animations;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Objects
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class BoxObject : public GeometryObject<N, Real_t> {
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIAS
    ////////////////////////////////////////////////////////////////////////////////
public:
    BoxObject() = delete;
    BoxObject(const JParams& jParams) { parseParameters(jParams); }

    virtual String   name() override { return String("BoxObject"); }
    virtual Real_t signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;

    virtual VecN getAABBMin() const override;
    virtual VecN getAABBMax() const override;

    void        setOriginalBox(const VecN& bMin, const VecN& bMax) { m_BoxMin = bMin; m_BoxMax = bMax; }
    const auto& originalBoxMin() const { return m_BoxMin; }
    const auto& originalBoxMax() const { return m_BoxMax; }
    auto getTransformedBoxMin() const { return this->transform(m_BoxMin); }
    auto getTransformedBoxMax() const { return this->transform(m_BoxMax); }

protected:
    virtual void parseParameters(const JParams& jParams) override;

    VecN m_BoxMin = VecN(-1.0);
    VecN m_BoxMax = VecN(1.0);
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class SphereObject : public GeometryObject<N, Real_t> {
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIAS
    ////////////////////////////////////////////////////////////////////////////////
public:
    SphereObject() = delete;
    SphereObject(const JParams& jParams) { this->parseParameters(jParams); }
    virtual String   name() override { if(N == 2) { return String("CircleObject"); } else { return String("SphereObject"); } }
    virtual Real_t signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class TorusObject : public GeometryObject<N, Real_t> {
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIAS
    ////////////////////////////////////////////////////////////////////////////////
public:
    TorusObject() = delete;
    TorusObject(const JParams& jParams) { parseParameters(jParams); }
    virtual String   name() override { return String("TorusObject"); }
    virtual Real_t signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
    void setRingRadius(Real_t ringRadius) { __NT_REQUIRE(ringRadius > 0); m_RingRadius = ringRadius; }
protected:
    virtual void parseParameters(const JParams& jParams) override;
    Real_t m_OuterRadius = Real_t(1);
    Real_t m_RingRadius  = Real_t(0.25);
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class Torus28Object : public TorusObject<N, Real_t> {
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIAS
    ////////////////////////////////////////////////////////////////////////////////
public:
    Torus28Object() = delete;
    Torus28Object(const JParams& jParams) : TorusObject<N, Real_t>(jParams) { this->parseParameters(jParams); }
    virtual Real_t signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class Torus2InfObject : public TorusObject<N, Real_t> {
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIAS
    ////////////////////////////////////////////////////////////////////////////////
public:
    Torus2InfObject() = delete;
    Torus2InfObject(const JParams& jParams) : TorusObject<N, Real_t>(jParams) { this->parseParameters(jParams); }
    virtual Real_t signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class Torus88Object : public TorusObject<N, Real_t> {
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIAS
    ////////////////////////////////////////////////////////////////////////////////
public:
    Torus88Object() = delete;
    Torus88Object(const JParams& jParams) : TorusObject<N, Real_t>(jParams) { this->parseParameters(jParams); }
    virtual Real_t signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class TorusInfInfObject : public TorusObject<N, Real_t> {
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIAS
    ////////////////////////////////////////////////////////////////////////////////
public:
    TorusInfInfObject() = delete;
    TorusInfInfObject(const JParams& jParams) : TorusObject<N, Real_t>(jParams) { this->parseParameters(jParams); }
    virtual Real_t signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class CylinderObject : public GeometryObject<N, Real_t> {
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIAS
    ////////////////////////////////////////////////////////////////////////////////
public:
    CylinderObject() = delete;
    CylinderObject(const JParams& jParams) { parseParameters(jParams); }
    virtual String   name() override { return String("CylinderObject"); }
    virtual Real_t signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
    void setRadius(Real_t radius) { m_Radius = radius; }
protected:
    virtual void parseParameters(const JParams& jParams) override;
    Real_t m_Radius = Real_t(0.5);
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class ConeObject : public GeometryObject<N, Real_t> {
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIAS
    ////////////////////////////////////////////////////////////////////////////////
public:
    ConeObject() = delete;
    ConeObject(const JParams& jParams) { parseParameters(jParams); }
    virtual String   name() override { return String("ConeObject"); }
    virtual Real_t signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
    void setBaseRadius(Real_t radius) { m_Radius = radius; }
protected:
    virtual void parseParameters(const JParams& jParams) override;
    Real_t m_Radius = Real_t(1);
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class PlaneObject : public GeometryObject<N, Real_t> {
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIAS
    ////////////////////////////////////////////////////////////////////////////////
public:
    PlaneObject() = delete;
    PlaneObject(const JParams& jParams) { parseParameters(jParams); }
    virtual String   name() override { return String("PlaneObject"); }
    virtual Real_t signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;

    void setNormal(const VecN& normal) { m_Normal = normal; }
    void setOffset(Real_t offset) { m_Offset = offset; }
protected:
    virtual void parseParameters(const JParams& jParams) override;
    VecN   m_Normal = VecN(0);
    Real_t m_Offset = Real_t(0);
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class TriangleObject : public GeometryObject<N, Real_t> {
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIAS
    ////////////////////////////////////////////////////////////////////////////////
public:
    TriangleObject() = delete;
    TriangleObject(const JParams& jParams) { parseParameters(jParams); }
    virtual String   name() override { return String("TriangleObject"); }
    virtual Real_t signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
    template<class IndexType> void setVertex(IndexType idx, const VecN& vertex) { m_Vertices[idx] = vertex; }
protected:
    virtual void parseParameters(const JParams& jParams) override;
    VecN m_Vertices[3];
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class HexagonObject : public GeometryObject<N, Real_t> {
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIAS
    ////////////////////////////////////////////////////////////////////////////////
public:
    HexagonObject() = delete;
    HexagonObject(const JParams& jParams) { this->parseParameters(jParams); }
    virtual String   name() override { return String("HexagonObject"); }
    virtual Real_t signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class TriangularPrismObject : public GeometryObject<N, Real_t> {
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIAS
    ////////////////////////////////////////////////////////////////////////////////
public:
    TriangularPrismObject() = delete;
    TriangularPrismObject(const JParams& jParams) { parseParameters(jParams); }
    virtual String   name() override { return String("TriangularPrismObject"); }
    virtual Real_t signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
    void setWidth(Real_t ratio) { m_Width = ratio; }
protected:
    virtual void parseParameters(const JParams& jParams) override;
    Real_t m_Width = Real_t(0.5);
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class HexagonalPrismObject : public GeometryObject<N, Real_t> {
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIAS
    ////////////////////////////////////////////////////////////////////////////////
public:
    HexagonalPrismObject() = delete;
    HexagonalPrismObject(const JParams& jParams) { parseParameters(jParams); }
    virtual String   name() override { return String("HexagonalPrismObject"); }
    virtual Real_t signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
    void             setWidth(Real_t ratio) { m_Width = ratio; }
protected:
    virtual void parseParameters(const JParams& jParams) override;
    Real_t m_Width = Real_t(0.5);
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class CapsuleObject : public GeometryObject<N, Real_t> {
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIAS
    ////////////////////////////////////////////////////////////////////////////////
public:
    CapsuleObject() = delete;
    CapsuleObject(const JParams& jParams) { parseParameters(jParams); }
    virtual String   name() override { return String("CapsuleObject"); }
    virtual Real_t signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
    void setRadius(Real_t radius) { m_Radius = radius; }
protected:
    virtual void parseParameters(const JParams& jParams) override;
    Real_t m_Radius = Real_t(0.5);
    VecN   m_Start  = VecN(0);
    VecN   m_End    = VecN(0);
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class EllipsoidObject : public GeometryObject<N, Real_t> {
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIAS
    ////////////////////////////////////////////////////////////////////////////////
public:
    EllipsoidObject() = delete;
    EllipsoidObject(const JParams& jParams) { parseParameters(jParams); }
    virtual String   name() override { return String("SphereObject"); }
    virtual Real_t signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;
    void setRadiusRatio(const VecN& scale) { m_RadiusRatio = scale; }
protected:
    virtual void parseParameters(const JParams& jParams) override;
    VecN m_RadiusRatio = VecN(1.0);
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class TriMeshObject : public GeometryObject<N, Real_t> {
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIAS
    ////////////////////////////////////////////////////////////////////////////////
public:
    TriMeshObject() = delete;
    TriMeshObject(const JParams& jParams) { parseParameters(jParams); }
    virtual String   name() override { return String("TriMeshObject"); }
    virtual Real_t signedDistance(const VecN&, bool bNegativeInside = true) const override;

    auto& meshFile() { return m_TriMeshFile; }
    auto& sdfStep() { return m_Step; }
    void computeSDF();
protected:
    virtual void parseParameters(const JParams& jParams) override;
    bool             m_bSDFGenerated = false;
    String           m_TriMeshFile   = String("");
    Real_t           m_Step          = Real_t(1.0 / 256.0);
    Array<N, Real_t> m_SDFData;
    Grid<N, Real_t>  m_Grid3D;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
enum CSGOperations {
    Overwrite,
    Union,
    Subtraction,
    Intersection,
    BlendExp,
    BlendPoly
};

enum DomainDeformation {
    None,
    Twist,
    CheapBend
};

template<Int N, class Real_t>
class CSGObject : public GeometryObject<N, Real_t> {
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIAS
    ////////////////////////////////////////////////////////////////////////////////
public:
    struct CSGData {
        SharedPtr<GeometryObject<N, Real_t>> obj = nullptr;
        CSGOperations                        op  = Union;
    };
    CSGObject() = delete;
    CSGObject(const JParams& jParams) { this->parseParameters(jParams); }
    virtual String   name() override { return String("CSGObject"); }
    virtual Real_t signedDistance(const VecN& ppos0, bool bNegativeInside = true) const override;

    void addObject(const CSGData& obj) { m_Objects.push_back(obj); }
    void addObject(const SharedPtr<GeometryObject<N, Real_t>>& obj, CSGOperations op = Union) { addObject({ obj, op }); }
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
} // end namespace NTCodeBase
