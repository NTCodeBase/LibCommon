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

#include <LibCommon/Array/ArrayHelpers.h>
#include <LibCommon/Geometry/MeshLoader.h>
#include <LibCommon/Utils/NumberHelpers.h>
#include <LibCommon/Utils/JSONHelpers.h>
#include <LibCommon/Geometry/GeometryObjects.h>
#include <LibCommon/Geometry/GeometryHelpers.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace GeometryObjects
{
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Base class
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
VecX<N, RealType> GeometryObject<N, RealType>::getAABBMin() const
{
    return transform(VecN(0)) - VecN(m_UniformScale) * std::sqrt(glm::compAdd(VecN(1.0)));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
VecX<N, RealType> GeometryObject<N, RealType>::getAABBMax() const
{
    return transform(VecN(0)) + VecN(m_UniformScale) * std::sqrt(glm::compAdd(VecN(1.0)));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
VecX<N, RealType> GeometryObject<N, RealType>::gradSignedDistance(const VecN& ppos, bool bNegativeInside /*= true*/, RealType dx /*= RealType(1e-4)*/) const
{
    if constexpr(N == 2)
    {
        RealType v00 = signedDistance(Vec2<RealType>(ppos[0] - dx, ppos[1] - dx), bNegativeInside);
        RealType v01 = signedDistance(Vec2<RealType>(ppos[0] - dx, ppos[1] + dx), bNegativeInside);
        RealType v10 = signedDistance(Vec2<RealType>(ppos[0] + dx, ppos[1] - dx), bNegativeInside);
        RealType v11 = signedDistance(Vec2<RealType>(ppos[0] + dx, ppos[1] + dx), bNegativeInside);

        RealType ddx0 = v10 - v00;
        RealType ddx1 = v11 - v01;

        RealType ddy0 = v01 - v00;
        RealType ddy1 = v11 - v10;

        return (Vec2<RealType>(ddx0, ddy0) + Vec2<RealType>(ddx1, ddy1)) * RealType(0.5);
    } else {
        __NT_REQUIRE_MSG(N == 3, "Array dimension != 2,3");

        RealType v000 = signedDistance(Vec3<RealType>(ppos[0] - dx, ppos[1] - dx, ppos[2] - dx), bNegativeInside);
        RealType v001 = signedDistance(Vec3<RealType>(ppos[0] - dx, ppos[1] - dx, ppos[2] + dx), bNegativeInside);
        RealType v010 = signedDistance(Vec3<RealType>(ppos[0] - dx, ppos[1] + dx, ppos[2] - dx), bNegativeInside);
        RealType v011 = signedDistance(Vec3<RealType>(ppos[0] - dx, ppos[1] + dx, ppos[2] + dx), bNegativeInside);

        RealType v100 = signedDistance(Vec3<RealType>(ppos[0] + dx, ppos[1] - dx, ppos[2] - dx), bNegativeInside);
        RealType v101 = signedDistance(Vec3<RealType>(ppos[0] + dx, ppos[1] - dx, ppos[2] + dx), bNegativeInside);
        RealType v110 = signedDistance(Vec3<RealType>(ppos[0] + dx, ppos[1] + dx, ppos[2] - dx), bNegativeInside);
        RealType v111 = signedDistance(Vec3<RealType>(ppos[0] + dx, ppos[1] + dx, ppos[2] + dx), bNegativeInside);

        RealType ddx00 = v100 - v000;
        RealType ddx10 = v110 - v010;
        RealType ddx01 = v101 - v001;
        RealType ddx11 = v111 - v011;

        RealType ddy00 = v010 - v000;
        RealType ddy10 = v110 - v100;
        RealType ddy01 = v011 - v001;
        RealType ddy11 = v111 - v101;

        RealType ddz00 = v001 - v000;
        RealType ddz10 = v101 - v100;
        RealType ddz01 = v011 - v010;
        RealType ddz11 = v111 - v110;

        return (Vec3<RealType>(ddx00, ddy00, ddz00) +
                Vec3<RealType>(ddx01, ddy01, ddz01) +
                Vec3<RealType>(ddx10, ddy10, ddz10) +
                Vec3<RealType>(ddx11, ddy11, ddz11)) * RealType(0.25);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void GeometryObject<N, RealType>::setTranslation(const VecN& translation)
{
    m_IntrinsicTranslation = translation;
    updateIntrinsicTransformation();
    updateTransformation();
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void GeometryObject<N, RealType>::setRotation(const VecX<N + 1, RealType>& rotation)
{
    if(rotation[N] != 0) {
        m_IntrinsicRotation = rotation;
        updateIntrinsicTransformation();
        updateTransformation();
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void GeometryObject<N, RealType>::setUniformScale(const RealType scaleVal)
{
    m_UniformScale = scaleVal;
    updateIntrinsicTransformation();
    updateTransformation();
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void GeometryObject<N, RealType>::resetTransformation()
{
    m_bTransformed                  = false;
    m_IntrinsicTranslation          = VecN(0);
    m_IntrinsicRotation             = VecNp1(VecN(1), 0);
    m_UniformScale                  = RealType(1.0);
    m_IntrinsicTransformationMatrix = MatNp1xNp1(1.0);
    m_TransformationMatrix          = MatNp1xNp1(1.0);
    m_InvTransformationMatrix       = MatNp1xNp1(1.0);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
bool GeometryObject<N, RealType>::updateTransformation(UInt frame /*= 0*/, RealType frameFraction /*= RealType(0)*/)
{
    if(frame > 0 && m_Animations.size() == 0) {
        return false;
    }
    ////////////////////////////////////////////////////////////////////////////////
    m_AnimationTransformationMatrix = MatNp1xNp1(1.0);
    for(auto& animation : m_Animations) {
        if(animation.isActive(frame)) {
            m_AnimationTransformationMatrix = animation.getTransformationMatrix(frame, frameFraction) * m_AnimationTransformationMatrix;
        } else {
            m_AnimationTransformationMatrix = animation.getInactiveTransformationMatrix(frame) * m_AnimationTransformationMatrix;
        }
    }
    m_TransformationMatrix    = m_AnimationTransformationMatrix * m_IntrinsicTransformationMatrix;
    m_InvTransformationMatrix = glm::inverse(m_TransformationMatrix);
    m_bTransformed = true;
    ////////////////////////////////////////////////////////////////////////////////
    return true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void GeometryObject<N, RealType>::updateIntrinsicTransformation()
{
    m_IntrinsicTransformationMatrix = MatNp1xNp1(1);
    m_IntrinsicTransformationMatrix = glm::scale(m_IntrinsicTransformationMatrix, VecN(m_UniformScale));
    m_IntrinsicTransformationMatrix = glm::rotate(m_IntrinsicTransformationMatrix, m_IntrinsicRotation[N], VecN(m_IntrinsicRotation));
    m_IntrinsicTransformationMatrix = glm::translate(m_IntrinsicTransformationMatrix, m_IntrinsicTranslation);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
VecX<N, RealType> GeometryObject<N, RealType>::transformAnimation(const VecN& ppos) const
{
    if(!m_bTransformed || m_Animations.size() == 0) {
        return ppos;
    } else {
        return VecN(m_AnimationTransformationMatrix * VecX<N + 1, RealType>(ppos, 1.0));
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
VecX<N, RealType> GeometryObject<N, RealType>::transform(const VecN& ppos) const
{
    if(!m_bTransformed) {
        return ppos;
    } else {
        return VecN(m_TransformationMatrix * VecX<N + 1, RealType>(ppos, 1.0));
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
VecX<N, RealType> GeometryObject<N, RealType>::invTransform(const VecN& ppos) const
{
    if(!m_bTransformed) {
        return ppos;
    } else {
        return VecN(m_InvTransformationMatrix * VecX<N + 1, RealType>(ppos, 1.0));
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

template<Int N, class RealType>
void GeometryObject<N, RealType>::parseParameters(const JParams& jParams)
{
    if(VecN translation; JSONHelpers::readVector(jParams, translation, "Translation")) {
        setTranslation(translation);
    }
    if(VecN rotationEulerAngles;
       JSONHelpers::readVector(jParams, rotationEulerAngles, "RotationEulerAngles") ||
       JSONHelpers::readVector(jParams, rotationEulerAngles, "RotationEulerAngle")) {
        setRotation(MathHelpers::EulerToAxisAngle(rotationEulerAngles, false));
    } else if(VecX<N + 1, RealType> rotationAxisAngle; JSONHelpers::readVector(jParams, rotationAxisAngle, "RotationAxisAngle")) {
        rotationAxisAngle[N] = glm::radians(rotationAxisAngle[N]);
        setRotation(rotationAxisAngle);
    }
    if(RealType scale; JSONHelpers::readValue(jParams, scale, "Scale")) {
        setUniformScale(scale);
    }
    ////////////////////////////////////////////////////////////////////////////////
    // animation data
    if(jParams.find("Animation") != jParams.end()) {
        for(auto& jAnimation   : jParams["Animation"]) {
            auto& animationObj = getAnimation().emplace_back(RigidBodyAnimation<N, RealType> {});

            if(bool bPeriodic; JSONHelpers::readBool(jAnimation, bPeriodic, "Periodic")) {
                animationObj.setPeriodic(bPeriodic);
            }
            if(UInt startFrame, endFrame;
               JSONHelpers::readValue(jAnimation, startFrame, "StartFrame") | JSONHelpers::readValue(jAnimation, endFrame, "EndFrame")) {
                animationObj.setAnimationRange(startFrame, endFrame);
            }

            __NT_REQUIRE(jAnimation.find("KeyFrames") != jAnimation.end());
            for(auto& jKeyFrame : jAnimation["KeyFrames"]) {
                auto& keyFrame = animationObj.addKeyFrame();
                __NT_REQUIRE(JSONHelpers::readValue(jKeyFrame, keyFrame.frame, "Frame"));

                // translation
                JSONHelpers::readVector(jKeyFrame, keyFrame.translation, "Translation");

                // rotation
                if(VecN rotationEulerAngles;
                   JSONHelpers::readVector(jKeyFrame, rotationEulerAngles, "RotationEulerAngles") ||
                   JSONHelpers::readVector(jKeyFrame, rotationEulerAngles, "RotationEulerAngle")) {
                    keyFrame.rotation = MathHelpers::EulerToAxisAngle(rotationEulerAngles, false, true);
                } else if(JSONHelpers::readVector(jKeyFrame, keyFrame.rotation, "RotationAxisAngle")) {
                    keyFrame.rotation = glm::radians(keyFrame.rotation);
                }
            }

            bool bCubicInterpolationTranslation = true;
            bool bCubicInterpolationRotation    = true;
            JSONHelpers::readBool(jAnimation, bCubicInterpolationTranslation, "CubicInterpolationTranslation");
            JSONHelpers::readBool(jAnimation, bCubicInterpolationRotation,    "CubicInterpolationRotation");
            animationObj.makeReady(bCubicInterpolationTranslation, bCubicInterpolationRotation);
        }
        updateTransformation();
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Specialized GeometryObjects
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
RealType BoxObject<N, RealType>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const
{
    auto     ppos = this->invTransform(ppos0);
    RealType mind = HugeReal();

    if(NumberHelpers::isInside(ppos, m_BoxMin, m_BoxMax)) {
        for(Int d = 0; d < N; ++d) {
            mind = MathHelpers::min(mind, ppos[d] - m_BoxMin[d], m_BoxMax[d] - ppos[d]);
        }
        mind = -mind;         // negative because inside
    } else {
        VecN cp;
        for(Int d = 0; d < N; ++d) {
            cp[d] = MathHelpers::max(MathHelpers::min(ppos[d], m_BoxMax[d]), m_BoxMin[d]);
        }
        mind = glm::length(ppos - cp);
    }

    if(!bNegativeInside) {
        mind = -mind;
    }

    return this->m_UniformScale * mind;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
VecX<N, RealType> BoxObject<N, RealType>::getAABBMin() const
{
    return this->transform(m_BoxMin) - VecN(this->m_UniformScale) * std::sqrt(glm::compAdd(VecN(1.0)));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
VecX<N, RealType> BoxObject<N, RealType>::getAABBMax() const
{
    return this->transform(m_BoxMax) + VecN(this->m_UniformScale) * std::sqrt(glm::compAdd(VecN(1.0)));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void BoxObject<N, RealType>::parseParameters(const JParams& jParams)
{
    GeometryObject<N, RealType>::parseParameters(jParams);
    if(VecN bMin, bMax; JSONHelpers::readVector(jParams, bMin, "BoxMin") && JSONHelpers::readVector(jParams, bMax, "BoxMax")) {
        setOriginalBox(bMin, bMax);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
RealType SphereObject<N, RealType>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const
{
    auto     ppos = invTransform(ppos0);
    RealType d    = this->m_UniformScale * (glm::length(ppos) - RealType(1.0));
    return bNegativeInside ? d : -d;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
RealType TorusObject<N, RealType>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const
{
    auto ppos = invTransform(ppos0);
    if constexpr(N == 2) {
        RealType q = std::abs(MathHelpers::norm2(ppos[0], ppos[1]) - m_OuterRadius);
        RealType d = this->m_UniformScale * (q - m_RingRadius);
        return bNegativeInside ? d : -d;
    } else {
        Vec2<RealType> q = Vec2<RealType>(MathHelpers::norm2(ppos[0], ppos[2]) - m_OuterRadius, ppos[1]);
        RealType       d = this->m_UniformScale * (glm::length(q) - m_RingRadius);
        return bNegativeInside ? d : -d;
    }
}

////////////////////////////////////////////////////////////////////////////////
template<Int N, class RealType>
void TorusObject<N, RealType>::parseParameters(const JParams& jParams)
{
    GeometryObject<N, RealType>::parseParameters(jParams);
    ////////////////////////////////////////////////////////////////////////////////
    RealType ringRadius;
    if(JSONHelpers::readValue(jParams, ringRadius, "RingRadius")) {
        setRingRadius(ringRadius);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
RealType Torus28Object<N, RealType>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const
{
    __NT_REQUIRE_MSG(N == 3, "Object dimension != 3");
    auto           ppos = invTransform(ppos0);
    Vec2<RealType> q    = Vec2<RealType>(MathHelpers::norm2(ppos[0], ppos[2]) - this->m_OuterRadius, ppos[1]);
    RealType       d    = this->m_UniformScale * (MathHelpers::norm8(q[0], q[1]) - this->m_RingRadius);
    return bNegativeInside ? d : -d;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
RealType Torus2InfObject<N, RealType>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const
{
    __NT_REQUIRE_MSG(N == 3, "Object dimension != 3");
    auto ppos = invTransform(ppos0);

    Vec2<RealType> q = Vec2<RealType>(MathHelpers::norm2(ppos[0], ppos[2]) - this->m_OuterRadius, ppos[1]);
    RealType       d = this->m_UniformScale * (MathHelpers::norm_inf(q[0], q[1]) - this->m_RingRadius);
    return bNegativeInside ? d : -d;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
RealType Torus88Object<N, RealType>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const
{
    auto ppos = invTransform(ppos0);
    if constexpr(N == 2) {
        RealType d = this->m_UniformScale * (std::abs(MathHelpers::norm8(ppos[0], ppos[1]) - this->m_OuterRadius) - this->m_RingRadius);
        return bNegativeInside ? d : -d;
    } else {
        Vec2<RealType> q = Vec2<RealType>(MathHelpers::norm8(ppos[0], ppos[2]) - this->m_OuterRadius, ppos[1]);
        RealType       d = this->m_UniformScale * (MathHelpers::norm8(q[0], q[1]) - this->m_RingRadius);
        return bNegativeInside ? d : -d;
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
RealType TorusInfInfObject<N, RealType>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const
{
    auto ppos = invTransform(ppos0);
    if constexpr(N == 2) {
        RealType d = this->m_UniformScale *
                     (std::abs(MathHelpers::norm_inf(ppos[0], ppos[1]) - this->m_OuterRadius) - this->m_RingRadius);
        return bNegativeInside ? d : -d;
    } else {
        Vec2<RealType> q = Vec2<RealType>(MathHelpers::norm_inf(ppos[0], ppos[2]) - this->m_OuterRadius, ppos[1]);
        RealType       d = this->m_UniformScale * (MathHelpers::norm_inf(q[0], q[1]) - this->m_RingRadius);
        return bNegativeInside ? d : -d;
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
RealType CylinderObject<N, RealType>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const
{
    __NT_REQUIRE_MSG(N == 3, "Object dimension != 3");
    auto     ppos = invTransform(ppos0);
    RealType d    = this->m_UniformScale * MathHelpers::max(MathHelpers::norm2(ppos[0], ppos[2]) - m_Radius, std::abs(ppos[1]) - RealType(1.0));
    return bNegativeInside ? d : -d;
}

////////////////////////////////////////////////////////////////////////////////
template<Int N, class RealType>
void CylinderObject<N, RealType>::parseParameters(const JParams& jParams)
{
    GeometryObject<N, RealType>::parseParameters(jParams);
    ////////////////////////////////////////////////////////////////////////////////
    RealType radius;
    if(JSONHelpers::readValue(jParams, radius, "Radius")) {
        setRadius(radius);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
RealType ConeObject<N, RealType>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const
{
    __NT_REQUIRE_MSG(N == 3, "Object dimension != 3");
    auto     ppos  = invTransform(ppos0);
    RealType theta = std::atan(m_Radius);     // radius / h, where h = 1
    RealType d1    = MathHelpers::norm2(ppos[0], ppos[2]) * cos(theta) - std::abs(RealType(1) - ppos[1]) * sin(theta);
    auto     d     = this->m_UniformScale * MathHelpers::max(d1, ppos[1] - RealType(1), -ppos[1]);
    return bNegativeInside ? d : -d;
}

////////////////////////////////////////////////////////////////////////////////
template<Int N, class RealType>
void ConeObject<N, RealType>::parseParameters(const JParams& jParams)
{
    GeometryObject<N, RealType>::parseParameters(jParams);
    ////////////////////////////////////////////////////////////////////////////////
    RealType baseRadius;
    if(JSONHelpers::readValue(jParams, baseRadius, "BaseRadius")) {
        setBaseRadius(baseRadius);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
RealType PlaneObject<N, RealType>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const
{
    auto     ppos = invTransform(ppos0);
    RealType d    = glm::dot(ppos, m_Normal) - m_Offset;
    return bNegativeInside ? d : -d;
}

////////////////////////////////////////////////////////////////////////////////
template<Int N, class RealType>
void PlaneObject<N, RealType>::parseParameters(const JParams& jParams)
{
    GeometryObject<N, RealType>::parseParameters(jParams);
    ////////////////////////////////////////////////////////////////////////////////
    VecN     normal;
    RealType offset;
    if(JSONHelpers::readVector(jParams, normal, "Normal")) {
        setNormal(normal);
    }
    if(JSONHelpers::readValue(jParams, offset, "Offset")) {
        setOffset(offset);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
RealType TriangleObject<N, RealType>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const
{
    __NT_REQUIRE_MSG(N == 2, "Object dimension != 2");
    auto ppos = invTransform(ppos0);
    if constexpr(N == 2) {
        auto p = VecX<N + 1, RealType>(ppos, 0);
        auto a = VecX<N + 1, RealType>(m_Vertices[0], 0);
        auto b = VecX<N + 1, RealType>(m_Vertices[1], 0);
        auto c = VecX<N + 1, RealType>(m_Vertices[2], 0);

        VecX<N + 1, RealType> ba = b - a; VecX<N + 1, RealType> pa = p - a;
        VecX<N + 1, RealType> cb = c - b; VecX<N + 1, RealType> pb = p - b;
        VecX<N + 1, RealType> ac = a - c; VecX<N + 1, RealType> pc = p - c;
        auto                  nor = glm::cross(ba, ac);

        auto sgn = [](auto val) -> int
                   {
                       return (val > 0) - (val < 0);
                   };

        if(sgn(glm::dot(glm::cross(ba, nor), pa)) +
           sgn(glm::dot(glm::cross(cb, nor), pb)) +
           sgn(glm::dot(glm::cross(ac, nor), pc)) < 2) {
            auto d = sqrt(std::min(std::min(glm::length2(ba * MathHelpers::clamp(dot(ba, pa) / glm::length2(ba), RealType(0), RealType(1.0)) - pa),
                                            glm::length2(cb * MathHelpers::clamp(dot(cb, pb) / glm::length2(cb), RealType(0), RealType(1.0)) - pb)),
                                   glm::length2(ac * MathHelpers::clamp(dot(ac, pc) / glm::length2(ac), RealType(0), RealType(1.0)) - pc)));
            return bNegativeInside ? d : -d;
        } else {
            auto d = -sqrt(std::min(std::min(glm::length2(ba * MathHelpers::clamp(dot(ba, pa) / glm::length2(ba), RealType(0), RealType(1.0)) - pa),
                                             glm::length2(cb * MathHelpers::clamp(dot(cb, pb) / glm::length2(cb), RealType(0), RealType(1.0)) - pb)),
                                    glm::length2(ac * MathHelpers::clamp(dot(ac, pc) / glm::length2(ac), RealType(0), RealType(1.0)) - pc)));
            return bNegativeInside ? d : -d;
        }
    } else {
        __NT_CALLED_TO_WRONG_PLACE
        return 0;
    }
}

////////////////////////////////////////////////////////////////////////////////
template<Int N, class RealType>
void TriangleObject<N, RealType>::parseParameters(const JParams& jParams)
{
    GeometryObject<N, RealType>::parseParameters(jParams);
    ////////////////////////////////////////////////////////////////////////////////
    VecN vertices[3];
    __NT_REQUIRE(JSONHelpers::readVector(jParams, vertices[0], "V0"));
    __NT_REQUIRE(JSONHelpers::readVector(jParams, vertices[1], "V1"));
    __NT_REQUIRE(JSONHelpers::readVector(jParams, vertices[2], "V2"));
    for(UInt i = 0; i < 3; ++i) {
        setVertex(i, vertices[i]);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
RealType HexagonObject<N, RealType>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const
{
    __NT_REQUIRE_MSG(N == 2, "Object dimension != 2");
    auto     ppos = invTransform(ppos0);
    RealType dx   = fabs(ppos[0]);
    RealType dy   = fabs(ppos[1]);
    RealType d    = this->m_UniformScale * (MathHelpers::max((dx * RealType(0.866025) + dy * RealType(0.5)), dy) - RealType(1.0));
    return bNegativeInside ? d : -d;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
RealType TriangularPrismObject<N, RealType>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const
{
    __NT_REQUIRE_MSG(N == 3, "Object dimension != 3");
    auto ppos = invTransform(ppos0);
    VecN q;
    for(Int d = 0; d < N; ++d) {
        q[d] = std::abs(ppos[d]);
    }
    RealType d = this->m_UniformScale * MathHelpers::max(q[0] - m_Width,
                                                         MathHelpers::max(q[2] * RealType(0.866025) + ppos[1] * RealType(0.5),
                                                                          -ppos[1]) - RealType(0.5));
    return bNegativeInside ? d : -d;
}

////////////////////////////////////////////////////////////////////////////////
template<Int N, class RealType>
void TriangularPrismObject<N, RealType>::parseParameters(const JParams& jParams)
{
    GeometryObject<N, RealType>::parseParameters(jParams);
    ////////////////////////////////////////////////////////////////////////////////
    RealType width;
    if(JSONHelpers::readValue(jParams, width, "Width")) {
        setWidth(width);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
RealType HexagonalPrismObject<N, RealType>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const
{
    __NT_REQUIRE_MSG(N == 3, "Object dimension != 3");
    auto ppos = invTransform(ppos0);
    VecN q;
    for(Int d = 0; d < N; ++d) {
        q[d] = std::abs(ppos[d]);
    }
    RealType d = this->m_UniformScale * MathHelpers::max(q[0] - m_Width,
                                                         MathHelpers::max(q[2] * RealType(0.866025) + q[1] * RealType(0.5), q[1]) - RealType(1.0));
    return bNegativeInside ? d : -d;
}

////////////////////////////////////////////////////////////////////////////////
template<Int N, class RealType>
void HexagonalPrismObject<N, RealType>::parseParameters(const JParams& jParams)
{
    GeometryObject<N, RealType>::parseParameters(jParams);
    ////////////////////////////////////////////////////////////////////////////////
    RealType width;
    if(JSONHelpers::readValue(jParams, width, "Width")) {
        setWidth(width);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
RealType CapsuleObject<N, RealType>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const
{
    auto ppos = invTransform(ppos0);
    VecN pa   = ppos - m_Start;
    VecN ba   = m_End - m_Start;

    RealType h = MathHelpers::clamp(glm::dot(pa, ba) / glm::dot(ba, ba), RealType(0.0), RealType(1.0));
    RealType d = this->m_UniformScale * (glm::length(pa - ba * h) - m_Radius);
    return bNegativeInside ? d : -d;
}

////////////////////////////////////////////////////////////////////////////////
template<Int N, class RealType>
void CapsuleObject<N, RealType>::parseParameters(const JParams& jParams)
{
    GeometryObject<N, RealType>::parseParameters(jParams);
    ////////////////////////////////////////////////////////////////////////////////
    if(RealType radius; JSONHelpers::readValue(jParams, radius, "Radius")) {
        setRadius(radius);
    }
    ////////////////////////////////////////////////////////////////////////////////
    m_Start[0] = RealType(-1);
    m_End[0]   = RealType(1);
    JSONHelpers::readValue(jParams, m_Start[0], "Start");
    JSONHelpers::readValue(jParams, m_End[0],   "End");
    m_Start[0] += m_Radius;
    m_End[0]   -= m_Radius;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
RealType EllipsoidObject<N, RealType>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const
{
    auto ppos = invTransform(ppos0);
    for(Int d = 0; d < N; ++d) {
        ppos[d] /= m_RadiusRatio[d];
    }
    RealType d = this->m_UniformScale * glm::compMin(m_RadiusRatio) * (glm::length(ppos) - RealType(1.0));
    return bNegativeInside ? d : -d;
}

////////////////////////////////////////////////////////////////////////////////
template<Int N, class RealType>
void EllipsoidObject<N, RealType>::parseParameters(const JParams& jParams)
{
    GeometryObject<N, RealType>::parseParameters(jParams);
    ////////////////////////////////////////////////////////////////////////////////
    VecN radiusRatio;
    if(JSONHelpers::readVector(jParams, radiusRatio, "RadiusRatio")) {
        setRadiusRatio(radiusRatio);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// sign distance field for triangle mesh
template<class RealType>
void computeSDFMesh(const StdVT<Vec3ui>& faces, const StdVT_Vec3<RealType>& vertices, const Vec3<RealType>& origin, RealType cellSize,
                    UInt ni, UInt nj, UInt nk, Array<3, RealType>& SDF, Int exactBand = 1)
{
    __NT_REQUIRE(ni > 0 && nj > 0 && nk > 0);

    SDF.resize(ni, nj, nk);
    SDF.assign(RealType(ni + nj + nk) * cellSize);                   // upper bound on distance
    Array3ui closest_tri(ni, nj, nk, 0xffffffff);

    // intersection_count(i,j,k) is # of tri intersections in (i-1,i]x{j}x{k}
    // we begin by initializing distances near the mesh, and figuring out intersection counts
    Array3ui intersectionCount(ni, nj, nk, 0u);

    for(UInt face = 0, faceEnd = static_cast<UInt>(faces.size()); face < faceEnd; ++face) {
        UInt p = faces[face][0];
        UInt q = faces[face][1];
        UInt r = faces[face][2];

        // coordinates in grid to high precision
        Vec3<RealType> fp = (vertices[p] - origin) / cellSize;
        Vec3<RealType> fq = (vertices[q] - origin) / cellSize;
        Vec3<RealType> fr = (vertices[r] - origin) / cellSize;

        // do distances nearby
        Int i0 = MathHelpers::clamp(static_cast<Int>(MathHelpers::min(fp[0], fq[0], fr[0])) - exactBand, 0, static_cast<Int>(ni - 1));
        Int i1 = MathHelpers::clamp(static_cast<Int>(MathHelpers::max(fp[0], fq[0], fr[0])) + exactBand + 1, 0, static_cast<Int>(ni - 1));
        Int j0 = MathHelpers::clamp(static_cast<Int>(MathHelpers::min(fp[1], fq[1], fr[1])) - exactBand, 0, static_cast<Int>(nj - 1));
        Int j1 = MathHelpers::clamp(static_cast<Int>(MathHelpers::max(fp[1], fq[1], fr[1])) + exactBand + 1, 0, static_cast<Int>(nj - 1));
        Int k0 = MathHelpers::clamp(static_cast<Int>(MathHelpers::min(fp[2], fq[2], fr[2])) - exactBand, 0, static_cast<Int>(nk - 1));
        Int k1 = MathHelpers::clamp(static_cast<Int>(MathHelpers::max(fp[2], fq[2], fr[2])) + exactBand + 1, 0, static_cast<Int>(nk - 1));

        Scheduler::parallel_for<Int>(i0, i1 + 1, j0, j1 + 1, k0, k1 + 1,
                                     [&](Int i, Int j, Int k)
                                     {
                                         Vec3<RealType> gx = Vec3<RealType>(i, j, k) * cellSize + origin;
                                         RealType d        = GeometryHelpers::point_triangle_distance(gx, vertices[p], vertices[q], vertices[r]);

                                         if(d < SDF(i, j, k)) {
                                             SDF(i, j, k)         = d;
                                             closest_tri(i, j, k) = face;
                                         }
                                     });

        // and do intersection counts
        j0 = MathHelpers::clamp(static_cast<Int>(std::ceil(MathHelpers::min(fp[1], fq[1], fr[1]))) - 10, 0, static_cast<Int>(nj - 1));
        j1 = MathHelpers::clamp(static_cast<Int>(std::floor(MathHelpers::max(fp[1], fq[1], fr[1]))) + 10, 0, static_cast<Int>(nj - 1));
        k0 = MathHelpers::clamp(static_cast<Int>(std::ceil(MathHelpers::min(fp[2], fq[2], fr[2]))) - 10, 0, static_cast<Int>(nk - 1));
        k1 = MathHelpers::clamp(static_cast<Int>(std::floor(MathHelpers::max(fp[2], fq[2], fr[2]))) + 10, 0, static_cast<Int>(nk - 1));

        for(Int k = k0; k <= k1; ++k) {
            for(Int j = j0; j <= j1; ++j) {
                RealType a, b, c;

                if(GeometryHelpers::point_in_triangle_2d(static_cast<RealType>(j), static_cast<RealType>(k), fp[1], fp[2], fq[1], fq[2], fr[1], fr[2], a, b, c)) {
                    // intersection i coordinate
                    RealType fi = a * fp[0] + b * fq[0] + c * fr[0];

                    // intersection is in (i_interval-1,i_interval]
                    Int i_interval = MathHelpers::max(static_cast<Int>(std::ceil(fi)), 0);

                    // we enlarge the first interval to include everything to the -x direction
                    // we ignore intersections that are beyond the +x side of the grid
                    if(i_interval < static_cast<Int>(ni)) {
                        ++intersectionCount(i_interval, j, k);
                    }
                }
            }
        }
    }             // end loop face

    // and now we fill in the rest of the distances with fast sweeping
    for(UInt pass = 0; pass < 2; ++pass) {
        GeometryHelpers::sweep(faces, vertices, SDF, closest_tri, origin, cellSize,  +1,  +1,  +1);
        GeometryHelpers::sweep(faces, vertices, SDF, closest_tri, origin, cellSize, -1,  -1,  -1);
        GeometryHelpers::sweep(faces, vertices, SDF, closest_tri, origin, cellSize,  +1,  +1, -1);
        GeometryHelpers::sweep(faces, vertices, SDF, closest_tri, origin, cellSize, -1,  -1,   +1);
        GeometryHelpers::sweep(faces, vertices, SDF, closest_tri, origin, cellSize,  +1, -1,   +1);
        GeometryHelpers::sweep(faces, vertices, SDF, closest_tri, origin, cellSize, -1,   +1, -1);
        GeometryHelpers::sweep(faces, vertices, SDF, closest_tri, origin, cellSize,  +1, -1,  -1);
        GeometryHelpers::sweep(faces, vertices, SDF, closest_tri, origin, cellSize, -1,   +1,  +1);
    }

    // then figure out signs (inside/outside) from intersection counts
    Scheduler::parallel_for<UInt>(0, nk,
                                  [&](UInt k)
                                  {
                                      for(UInt j = 0; j < nj; ++j) {
                                          UInt total_count = 0;

                                          for(UInt i = 0; i < ni; ++i) {
                                              total_count += intersectionCount(i, j, k);

                                              if(total_count & 1) {                 // if parity of intersections so far is odd,
                                                  SDF(i, j, k) = -SDF(i, j, k);     // we are inside the mesh
                                              }
                                          }
                                      }
                                  });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
RealType TriMeshObject<N, RealType>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const
{
    if constexpr(N == 2) {
        __NT_UNUSED(ppos0);
        __NT_UNUSED(bNegativeInside);
        return 0;
    } else {
        __NT_REQUIRE(m_bSDFGenerated);
        auto     ppos    = invTransform(ppos0);
        auto     gridPos = m_Grid3D.getGridCoordinate(ppos);
        RealType d       = this->m_UniformScale * ArrayHelpers::interpolateValueLinear(gridPos, m_SDFData);
        return bNegativeInside ? d : -d;
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void TriMeshObject<N, RealType>::computeSDF()
{
    if constexpr(N == 3) {
        ////////////////////////////////////////////////////////////////////////////////
        // Load mesh
        MeshLoader meshLoader;
        __NT_REQUIRE(meshLoader.loadMesh(m_TriMeshFile));
        meshLoader.scaleToBox();

        ////////////////////////////////////////////////////////////////////////////////
        // Init grid and array of data
        m_Grid3D.setGrid(VecN(meshLoader.getAABBMin()) - VecN(RealType(3.0) * m_Step),
                         VecN(meshLoader.getAABBMax()) + VecN(RealType(3.0) * m_Step),
                         m_Step);

        StdVT_Vec3<RealType> vertexList(meshLoader.getNVertices());
        StdVT<Vec3ui>        faceList(meshLoader.getNFaces());

        std::memcpy(vertexList.data(), meshLoader.getVertices().data(), meshLoader.getVertices().size() * sizeof(RealType));
        std::memcpy(faceList.data(),   meshLoader.getFaces().data(),    meshLoader.getFaces().size() * sizeof(UInt));

        ////////////////////////////////////////////////////////////////////////////////
        // Compute SDF data
        computeSDFMesh(faceList, vertexList, VecN(meshLoader.getAABBMin()),
                       m_Step, m_Grid3D.getNCells()[0], m_Grid3D.getNCells()[1], m_Grid3D.getNCells()[2], m_SDFData);
        m_bSDFGenerated = true;
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void TriMeshObject<N, RealType>::parseParameters(const JParams& jParams)
{
    GeometryObject<N, RealType>::parseParameters(jParams);
    ////////////////////////////////////////////////////////////////////////////////
    __NT_REQUIRE(JSONHelpers::readValue(jParams, meshFile(), "MeshFile"));
    JSONHelpers::readValue(jParams, sdfStep(), "SDFStep");
    computeSDF();
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
RealType CSGObject<N, RealType>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const
{
    if(m_Objects.size() == 0) {
        return HugeReal();
    }

    auto     ppos = domainDeform(transform(ppos0));
    RealType sd   = m_Objects[0].obj->signedDistance(ppos, bNegativeInside);

    for(size_t i = 1; i < m_Objects.size(); ++i) {
        auto& csgObj = m_Objects[i];
        switch(csgObj.op) {
            case Overwrite:
                sd = csgObj.obj->signedDistance(ppos);
                break;
            case Union:
                sd = MathHelpers::min(sd, csgObj.obj->signedDistance(ppos, bNegativeInside));
                break;
            case Subtraction:
                sd = MathHelpers::max(sd, -csgObj.obj->signedDistance(ppos, bNegativeInside));
                break;
            case Intersection:
                sd = MathHelpers::max(sd, csgObj.obj->signedDistance(ppos, bNegativeInside));
                break;
            case BlendExp:
                sd = MathHelpers::smin_exp(sd, csgObj.obj->signedDistance(ppos, bNegativeInside));
                break;
            case BlendPoly:
                sd = MathHelpers::smin_poly(sd, csgObj.obj->signedDistance(ppos, bNegativeInside));
                break;
            default:
                ;
        }
    }

    return sd;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
VecX<N, RealType> CSGObject<N, RealType>::domainDeform(const VecN& ppos) const
{
    switch(m_DeformOp) {
        case None:
            return ppos;
        case Twist:
            return twist(ppos);
        case CheapBend:
            return cheapBend(ppos);
        default:
            return ppos;
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
VecX<N, RealType> CSGObject<N, RealType>::twist(const VecN& ppos) const
{
    __NT_UNUSED(ppos);
    //RealType         c = cos(RealType(5.0) * ppos.z);
    //RealType         s = sin(RealType(5.0) * ppos.z);
    //Mat2x2<RealType> m = Mat2x2<RealType>(c, -s, s, c);
    //Mat2x2<RealType> m = Mat2x2<RealType>(c, s, -s, c);

    //return VecN(m * VecN(ppos[0], ppos[1]), ppos.z);
    return VecN(0);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
VecX<N, RealType> CSGObject<N, RealType>::cheapBend(const VecN& ppos) const
{
    RealType c = cos(RealType(0.5) * ppos.y);
    RealType s = sin(RealType(0.5) * ppos.y);
    //Mat2x2<RealType> m = Mat2x2<RealType>(c, -s, s, c);
    Mat2x2<RealType> m = Mat2x2<RealType>(c, s, -s, c);

    //return VecN(m * VecN(ppos[0], ppos[1]), ppos.z);
    return VecN(0);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(GeometryObject)
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(BoxObject)
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(SphereObject)
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(TorusObject)
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(Torus28Object)
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(Torus2InfObject)
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(Torus88Object)
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(TorusInfInfObject)
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(CylinderObject)
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(ConeObject)
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(PlaneObject)
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(TriangleObject)
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(HexagonObject)
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(TriangularPrismObject)
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(HexagonalPrismObject)
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(CapsuleObject)
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(EllipsoidObject)
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(TriMeshObject)
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(CSGObject)

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
}   // end namespace GeometryObjects
