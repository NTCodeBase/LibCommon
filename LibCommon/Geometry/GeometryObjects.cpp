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

#include <LibCommon/Array/ArrayHelpers.h>
#include <LibCommon/Geometry/MeshLoader.h>
#include <LibCommon/Utils/NumberHelpers.h>
#include <LibCommon/Utils/JSONHelpers.h>
#include <LibCommon/Geometry/GeometryObjects.h>
#include <LibCommon/Geometry/GeometryHelpers.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NTCodeBase {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Base class
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
VecX<N, Real_t> GeometryObject<N, Real_t>::getAABBMin() const {
    return transform(VecN(0)) - VecN(m_UniformScale) * std::sqrt(glm::compAdd(VecN(1.0)));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
VecX<N, Real_t> GeometryObject<N, Real_t>::getAABBMax() const {
    return transform(VecN(0)) + VecN(m_UniformScale) * std::sqrt(glm::compAdd(VecN(1.0)));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
VecX<N, Real_t> GeometryObject<N, Real_t>::gradSignedDistance(const VecN& ppos, bool bNegativeInside /*= true*/, Real_t dx /*= Real_t(1e-4)*/) const {
    if constexpr(N == 2)
    {
        Real_t v00 = signedDistance(Vec2<Real_t>(ppos[0] - dx, ppos[1] - dx), bNegativeInside);
        Real_t v01 = signedDistance(Vec2<Real_t>(ppos[0] - dx, ppos[1] + dx), bNegativeInside);
        Real_t v10 = signedDistance(Vec2<Real_t>(ppos[0] + dx, ppos[1] - dx), bNegativeInside);
        Real_t v11 = signedDistance(Vec2<Real_t>(ppos[0] + dx, ppos[1] + dx), bNegativeInside);

        Real_t ddx0 = v10 - v00;
        Real_t ddx1 = v11 - v01;

        Real_t ddy0 = v01 - v00;
        Real_t ddy1 = v11 - v10;

        return (Vec2<Real_t>(ddx0, ddy0) + Vec2<Real_t>(ddx1, ddy1)) * Real_t(0.5);
    } else {
        NT_REQUIRE_MSG(N == 3, "Array dimension != 2,3");

        Real_t v000 = signedDistance(Vec3<Real_t>(ppos[0] - dx, ppos[1] - dx, ppos[2] - dx), bNegativeInside);
        Real_t v001 = signedDistance(Vec3<Real_t>(ppos[0] - dx, ppos[1] - dx, ppos[2] + dx), bNegativeInside);
        Real_t v010 = signedDistance(Vec3<Real_t>(ppos[0] - dx, ppos[1] + dx, ppos[2] - dx), bNegativeInside);
        Real_t v011 = signedDistance(Vec3<Real_t>(ppos[0] - dx, ppos[1] + dx, ppos[2] + dx), bNegativeInside);

        Real_t v100 = signedDistance(Vec3<Real_t>(ppos[0] + dx, ppos[1] - dx, ppos[2] - dx), bNegativeInside);
        Real_t v101 = signedDistance(Vec3<Real_t>(ppos[0] + dx, ppos[1] - dx, ppos[2] + dx), bNegativeInside);
        Real_t v110 = signedDistance(Vec3<Real_t>(ppos[0] + dx, ppos[1] + dx, ppos[2] - dx), bNegativeInside);
        Real_t v111 = signedDistance(Vec3<Real_t>(ppos[0] + dx, ppos[1] + dx, ppos[2] + dx), bNegativeInside);

        Real_t ddx00 = v100 - v000;
        Real_t ddx10 = v110 - v010;
        Real_t ddx01 = v101 - v001;
        Real_t ddx11 = v111 - v011;

        Real_t ddy00 = v010 - v000;
        Real_t ddy10 = v110 - v100;
        Real_t ddy01 = v011 - v001;
        Real_t ddy11 = v111 - v101;

        Real_t ddz00 = v001 - v000;
        Real_t ddz10 = v101 - v100;
        Real_t ddz01 = v011 - v010;
        Real_t ddz11 = v111 - v110;

        return (Vec3<Real_t>(ddx00, ddy00, ddz00) +
                Vec3<Real_t>(ddx01, ddy01, ddz01) +
                Vec3<Real_t>(ddx10, ddy10, ddz10) +
                Vec3<Real_t>(ddx11, ddy11, ddz11)) * Real_t(0.25);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void GeometryObject<N, Real_t>::setTranslation(const VecN& translation) {
    m_IntrinsicTranslation = translation;
    updateIntrinsicTransformation();
    updateTransformation();
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void GeometryObject<N, Real_t>::setRotation(const VecNp1& rotation) {
    if(rotation[N] != 0) {
        m_IntrinsicRotation = rotation;
        updateIntrinsicTransformation();
        updateTransformation();
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void GeometryObject<N, Real_t>::setUniformScale(const Real_t scaleVal) {
    m_UniformScale = scaleVal;
    updateIntrinsicTransformation();
    updateTransformation();
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void GeometryObject<N, Real_t>::resetTransformation() {
    m_bTransformed                  = false;
    m_IntrinsicTranslation          = VecN(0);
    m_IntrinsicRotation             = VecNp1(VecN(1), 0);
    m_UniformScale                  = Real_t(1.0);
    m_IntrinsicTransformationMatrix = MatNp1xNp1(1.0);
    m_TransformationMatrix          = MatNp1xNp1(1.0);
    m_InvTransformationMatrix       = MatNp1xNp1(1.0);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
bool GeometryObject<N, Real_t>::updateTransformation(UInt frame /*= 0*/, Real_t frameFraction /*= Real_t(0)*/) {
    if(m_bDoneTransformation || (frame > 0 && m_Animations.size() == 0)) {
        return false;
    }
    ////////////////////////////////////////////////////////////////////////////////
    m_PrevAnimationTransformationMatrix = m_AnimationTransformationMatrix;
    m_AnimationTransformationMatrix     = MatNp1xNp1(1.0);
    bool bDone = (m_Animations.size() > 0);
    for(auto& animation : m_Animations) {
        if(animation.isActive(frame)) {
            m_AnimationTransformationMatrix = animation.getTransformationMatrix(frame, frameFraction) * m_AnimationTransformationMatrix;
        } else {
            m_AnimationTransformationMatrix = animation.getInactiveTransformationMatrix(frame) * m_AnimationTransformationMatrix;
        }
        if(!animation.doneAnimation(frame)) {
            bDone = false;
        }
    }
    m_InvAnimationTransformationMatrix = glm::inverse(m_AnimationTransformationMatrix);
    ////////////////////////////////////////////////////////////////////////////////
    m_PrevTransformationMatrix = m_TransformationMatrix;
    m_TransformationMatrix     = m_AnimationTransformationMatrix * m_IntrinsicTransformationMatrix;
    m_InvTransformationMatrix  = glm::inverse(m_TransformationMatrix);
    m_bTransformed        = true;
    m_bDoneTransformation = bDone;
    ////////////////////////////////////////////////////////////////////////////////
    return true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void GeometryObject<N, Real_t>::updateIntrinsicTransformation() {
    m_IntrinsicTransformationMatrix = MatNp1xNp1(1);
    m_IntrinsicTransformationMatrix = glm::translate(m_IntrinsicTransformationMatrix, m_IntrinsicTranslation);
    m_IntrinsicTransformationMatrix = glm::rotate(m_IntrinsicTransformationMatrix, m_IntrinsicRotation[N], VecN(m_IntrinsicRotation));
    m_IntrinsicTransformationMatrix = glm::scale(m_IntrinsicTransformationMatrix, VecN(m_UniformScale));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
VecX<N, Real_t> GeometryObject<N, Real_t>::transformAnimation(const VecN& ppos) const {
    if(!m_bTransformed || m_Animations.size() == 0) {
        return ppos;
    } else {
        return VecN(m_AnimationTransformationMatrix * VecNp1(ppos, 1.0));
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
VecX<N, Real_t> GeometryObject<N, Real_t>::invTransformAnimation(const VecN& ppos) const {
    if(!m_bTransformed || m_Animations.size() == 0) {
        return ppos;
    } else {
        return VecN(m_InvAnimationTransformationMatrix * VecNp1(ppos, 1.0));
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
VecX<N, Real_t> GeometryObject<N, Real_t>::transform(const VecN& ppos) const {
    if(!m_bTransformed) {
        return ppos;
    } else {
        return VecN(m_TransformationMatrix * VecNp1(ppos, 1.0));
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
VecX<N, Real_t> GeometryObject<N, Real_t>::invTransform(const VecN& ppos) const {
    if(!m_bTransformed) {
        return ppos;
    } else {
        return VecN(m_InvTransformationMatrix * VecNp1(ppos, 1.0));
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

template<Int N, class Real_t>
void GeometryObject<N, Real_t>::parseParameters(const JParams& jParams) {
    if(VecN translation; JSONHelpers::readVector(jParams, translation, "Translation")) {
        setTranslation(translation);
    }
    if(VecN rotationEulerAngles;
       JSONHelpers::readVector(jParams, rotationEulerAngles, "RotationEulerAngles") ||
       JSONHelpers::readVector(jParams, rotationEulerAngles, "RotationEulerAngle")) {
        setRotation(MathHelpers::EulerToAxisAngle(rotationEulerAngles, false));
    } else if(VecNp1 rotationAxisAngle; JSONHelpers::readVector(jParams, rotationAxisAngle, "RotationAxisAngle")) {
        rotationAxisAngle[N] = glm::radians(rotationAxisAngle[N]);
        setRotation(rotationAxisAngle);
    }
    if(Real_t scale; JSONHelpers::readValue(jParams, scale, "Scale")) {
        setUniformScale(scale);
    }
    ////////////////////////////////////////////////////////////////////////////////
    // animation data
    if(jParams.find("Animation") != jParams.end()) {
        for(auto& jAnimation   : jParams["Animation"]) {
            auto& animationObj = getAnimation().emplace_back(RigidBodyAnimation<N, Real_t> {});

            if(bool bPeriodic; JSONHelpers::readBool(jAnimation, bPeriodic, "Periodic")) {
                animationObj.periodic() = bPeriodic;
            }

            NT_REQUIRE(jAnimation.find("KeyFrames") != jAnimation.end());
            for(auto& jKeyFrame : jAnimation["KeyFrames"]) {
                auto& keyFrame = animationObj.addKeyFrame();
                NT_REQUIRE(JSONHelpers::readValue(jKeyFrame, keyFrame.frame, "Frame"));

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
template<Int N, class Real_t>
Real_t BoxObject<N, Real_t>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const {
    auto   ppos = this->invTransform(ppos0);
    Real_t mind = HugeReal();

    if(NumberHelpers::isInside(ppos, m_BoxMin, m_BoxMax)) {
        for(Int d = 0; d < N; ++d) {
            mind = MathHelpers::min(mind, ppos[d] - m_BoxMin[d], m_BoxMax[d] - ppos[d]);
        }
        mind = -mind; // negative because inside
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
template<Int N, class Real_t>
VecX<N, Real_t> BoxObject<N, Real_t>::getAABBMin() const {
    return this->transform(m_BoxMin) - VecN(this->m_UniformScale) * std::sqrt(glm::compAdd(VecN(1.0)));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
VecX<N, Real_t> BoxObject<N, Real_t>::getAABBMax() const {
    return this->transform(m_BoxMax) + VecN(this->m_UniformScale) * std::sqrt(glm::compAdd(VecN(1.0)));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void BoxObject<N, Real_t>::parseParameters(const JParams& jParams) {
    GeometryObject<N, Real_t>::parseParameters(jParams);
    if(VecN bMin, bMax; JSONHelpers::readVector(jParams, bMin, "BoxMin") && JSONHelpers::readVector(jParams, bMax, "BoxMax")) {
        setOriginalBox(bMin, bMax);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
Real_t SphereObject<N, Real_t>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const {
    auto   ppos = this->invTransform(ppos0);
    Real_t d    = this->m_UniformScale * (glm::length(ppos) - Real_t(1.0));
    return bNegativeInside ? d : -d;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
Real_t TorusObject<N, Real_t>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const {
    auto ppos = this->invTransform(ppos0);
    if constexpr(N == 2) {
        Real_t q = std::abs(MathHelpers::norm2(ppos[0], ppos[1]) - m_OuterRadius);
        Real_t d = this->m_UniformScale * (q - m_RingRadius);
        return bNegativeInside ? d : -d;
    } else {
        Vec2<Real_t> q = Vec2<Real_t>(MathHelpers::norm2(ppos[0], ppos[2]) - m_OuterRadius, ppos[1]);
        Real_t       d = this->m_UniformScale * (glm::length(q) - m_RingRadius);
        return bNegativeInside ? d : -d;
    }
}

////////////////////////////////////////////////////////////////////////////////
template<Int N, class Real_t>
void TorusObject<N, Real_t>::parseParameters(const JParams& jParams) {
    GeometryObject<N, Real_t>::parseParameters(jParams);
    ////////////////////////////////////////////////////////////////////////////////
    Real_t ringRadius;
    if(JSONHelpers::readValue(jParams, ringRadius, "RingRadius")) {
        setRingRadius(ringRadius);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
Real_t Torus28Object<N, Real_t>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const {
    NT_REQUIRE_MSG(N == 3, "Object dimension != 3");
    auto         ppos = this->invTransform(ppos0);
    Vec2<Real_t> q    = Vec2<Real_t>(MathHelpers::norm2(ppos[0], ppos[2]) - this->m_OuterRadius, ppos[1]);
    Real_t       d    = this->m_UniformScale * (MathHelpers::norm8(q[0], q[1]) - this->m_RingRadius);
    return bNegativeInside ? d : -d;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
Real_t Torus2InfObject<N, Real_t>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const {
    NT_REQUIRE_MSG(N == 3, "Object dimension != 3");
    auto ppos = this->invTransform(ppos0);

    Vec2<Real_t> q = Vec2<Real_t>(MathHelpers::norm2(ppos[0], ppos[2]) - this->m_OuterRadius, ppos[1]);
    Real_t       d = this->m_UniformScale * (MathHelpers::norm_inf(q[0], q[1]) - this->m_RingRadius);
    return bNegativeInside ? d : -d;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
Real_t Torus88Object<N, Real_t>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const {
    auto ppos = this->invTransform(ppos0);
    if constexpr(N == 2) {
        Real_t d = this->m_UniformScale * (std::abs(MathHelpers::norm8(ppos[0], ppos[1]) - this->m_OuterRadius) - this->m_RingRadius);
        return bNegativeInside ? d : -d;
    } else {
        Vec2<Real_t> q = Vec2<Real_t>(MathHelpers::norm8(ppos[0], ppos[2]) - this->m_OuterRadius, ppos[1]);
        Real_t       d = this->m_UniformScale * (MathHelpers::norm8(q[0], q[1]) - this->m_RingRadius);
        return bNegativeInside ? d : -d;
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
Real_t TorusInfInfObject<N, Real_t>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const {
    auto ppos = this->invTransform(ppos0);
    if constexpr(N == 2) {
        Real_t d = this->m_UniformScale *
                   (std::abs(MathHelpers::norm_inf(ppos[0], ppos[1]) - this->m_OuterRadius) - this->m_RingRadius);
        return bNegativeInside ? d : -d;
    } else {
        Vec2<Real_t> q = Vec2<Real_t>(MathHelpers::norm_inf(ppos[0], ppos[2]) - this->m_OuterRadius, ppos[1]);
        Real_t       d = this->m_UniformScale * (MathHelpers::norm_inf(q[0], q[1]) - this->m_RingRadius);
        return bNegativeInside ? d : -d;
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
Real_t CylinderObject<N, Real_t>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const {
    NT_REQUIRE_MSG(N == 3, "Object dimension != 3");
    auto   ppos = this->invTransform(ppos0);
    Real_t d    = this->m_UniformScale * MathHelpers::max(MathHelpers::norm2(ppos[0], ppos[2]) - m_Radius, std::abs(ppos[1]) - Real_t(1.0));
    return bNegativeInside ? d : -d;
}

////////////////////////////////////////////////////////////////////////////////
template<Int N, class Real_t>
void CylinderObject<N, Real_t>::parseParameters(const JParams& jParams) {
    GeometryObject<N, Real_t>::parseParameters(jParams);
    ////////////////////////////////////////////////////////////////////////////////
    Real_t radius;
    if(JSONHelpers::readValue(jParams, radius, "Radius")) {
        setRadius(radius);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
Real_t ConeObject<N, Real_t>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const {
    NT_REQUIRE_MSG(N == 3, "Object dimension != 3");
    auto   ppos  = this->invTransform(ppos0);
    Real_t theta = std::atan(m_Radius); // radius / h, where h = 1
    Real_t d1    = MathHelpers::norm2(ppos[0], ppos[2]) * cos(theta) - std::abs(Real_t(1) - ppos[1]) * sin(theta);
    auto   d     = this->m_UniformScale * MathHelpers::max(d1, ppos[1] - Real_t(1), -ppos[1]);
    return bNegativeInside ? d : -d;
}

////////////////////////////////////////////////////////////////////////////////
template<Int N, class Real_t>
void ConeObject<N, Real_t>::parseParameters(const JParams& jParams) {
    GeometryObject<N, Real_t>::parseParameters(jParams);
    ////////////////////////////////////////////////////////////////////////////////
    Real_t baseRadius;
    if(JSONHelpers::readValue(jParams, baseRadius, "BaseRadius")) {
        setBaseRadius(baseRadius);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
Real_t PlaneObject<N, Real_t>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const {
    auto   ppos = this->invTransform(ppos0);
    Real_t d    = glm::dot(ppos, m_Normal) - m_Offset;
    return bNegativeInside ? d : -d;
}

////////////////////////////////////////////////////////////////////////////////
template<Int N, class Real_t>
void PlaneObject<N, Real_t>::parseParameters(const JParams& jParams) {
    GeometryObject<N, Real_t>::parseParameters(jParams);
    ////////////////////////////////////////////////////////////////////////////////
    VecN   normal;
    Real_t offset;
    if(JSONHelpers::readVector(jParams, normal, "Normal")) {
        setNormal(normal);
    }
    if(JSONHelpers::readValue(jParams, offset, "Offset")) {
        setOffset(offset);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
Real_t TriangleObject<N, Real_t>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const {
    NT_REQUIRE_MSG(N == 2, "Object dimension != 2");
    auto ppos = this->invTransform(ppos0);
    if constexpr(N == 2) {
        auto p = VecNp1(ppos, 0);
        auto a = VecNp1(m_Vertices[0], 0);
        auto b = VecNp1(m_Vertices[1], 0);
        auto c = VecNp1(m_Vertices[2], 0);

        VecNp1 ba = b - a; VecNp1 pa = p - a;
        VecNp1 cb = c - b; VecNp1 pb = p - b;
        VecNp1 ac = a - c; VecNp1 pc = p - c;
        auto   nor = glm::cross(ba, ac);

        auto sgn = [](auto val) -> int {
                       return (val > 0) - (val < 0);
                   };

        if(sgn(glm::dot(glm::cross(ba, nor), pa)) +
           sgn(glm::dot(glm::cross(cb, nor), pb)) +
           sgn(glm::dot(glm::cross(ac, nor), pc)) < 2) {
            auto d = std::sqrt(std::min(std::min(glm::length2(ba * MathHelpers::clamp(dot(ba, pa) / glm::length2(ba), Real_t(0), Real_t(1.0)) - pa),
                                                 glm::length2(cb * MathHelpers::clamp(dot(cb, pb) / glm::length2(cb), Real_t(0), Real_t(1.0)) - pb)),
                                        glm::length2(ac * MathHelpers::clamp(dot(ac, pc) / glm::length2(ac), Real_t(0), Real_t(1.0)) - pc)));
            return bNegativeInside ? d : -d;
        } else {
            auto d = -sqrt(std::min(std::min(glm::length2(ba * MathHelpers::clamp(dot(ba, pa) / glm::length2(ba), Real_t(0), Real_t(1.0)) - pa),
                                             glm::length2(cb * MathHelpers::clamp(dot(cb, pb) / glm::length2(cb), Real_t(0), Real_t(1.0)) - pb)),
                                    glm::length2(ac * MathHelpers::clamp(dot(ac, pc) / glm::length2(ac), Real_t(0), Real_t(1.0)) - pc)));
            return bNegativeInside ? d : -d;
        }
    } else {
        NT_CALLED_TO_WRONG_PLACE
        return 0;
    }
}

////////////////////////////////////////////////////////////////////////////////
template<Int N, class Real_t>
void TriangleObject<N, Real_t>::parseParameters(const JParams& jParams) {
    GeometryObject<N, Real_t>::parseParameters(jParams);
    ////////////////////////////////////////////////////////////////////////////////
    VecN vertices[3];
    NT_REQUIRE(JSONHelpers::readVector(jParams, vertices[0], "V0"));
    NT_REQUIRE(JSONHelpers::readVector(jParams, vertices[1], "V1"));
    NT_REQUIRE(JSONHelpers::readVector(jParams, vertices[2], "V2"));
    for(UInt i = 0; i < 3; ++i) {
        setVertex(i, vertices[i]);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
Real_t HexagonObject<N, Real_t>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const {
    NT_REQUIRE_MSG(N == 2, "Object dimension != 2");
    auto   ppos = this->invTransform(ppos0);
    Real_t dx   = std::abs(ppos[0]);
    Real_t dy   = std::abs(ppos[1]);
    Real_t d    = this->m_UniformScale * (MathHelpers::max((dx * Real_t(0.866025) + dy * Real_t(0.5)), dy) - Real_t(1.0));
    return bNegativeInside ? d : -d;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
Real_t TriangularPrismObject<N, Real_t>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const {
    NT_REQUIRE_MSG(N == 3, "Object dimension != 3");
    auto ppos = this->invTransform(ppos0);
    VecN q;
    for(Int d = 0; d < N; ++d) {
        q[d] = std::abs(ppos[d]);
    }
    Real_t d = this->m_UniformScale * MathHelpers::max(q[0] - m_Width,
                                                       MathHelpers::max(q[2] * Real_t(0.866025) + ppos[1] * Real_t(0.5),
                                                                        -ppos[1]) - Real_t(0.5));
    return bNegativeInside ? d : -d;
}

////////////////////////////////////////////////////////////////////////////////
template<Int N, class Real_t>
void TriangularPrismObject<N, Real_t>::parseParameters(const JParams& jParams) {
    GeometryObject<N, Real_t>::parseParameters(jParams);
    ////////////////////////////////////////////////////////////////////////////////
    Real_t width;
    if(JSONHelpers::readValue(jParams, width, "Width")) {
        setWidth(width);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
Real_t HexagonalPrismObject<N, Real_t>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const {
    NT_REQUIRE_MSG(N == 3, "Object dimension != 3");
    auto ppos = this->invTransform(ppos0);
    VecN q;
    for(Int d = 0; d < N; ++d) {
        q[d] = std::abs(ppos[d]);
    }
    Real_t d = this->m_UniformScale * MathHelpers::max(q[0] - m_Width,
                                                       MathHelpers::max(q[2] * Real_t(0.866025) + q[1] * Real_t(0.5), q[1]) - Real_t(1.0));
    return bNegativeInside ? d : -d;
}

////////////////////////////////////////////////////////////////////////////////
template<Int N, class Real_t>
void HexagonalPrismObject<N, Real_t>::parseParameters(const JParams& jParams) {
    GeometryObject<N, Real_t>::parseParameters(jParams);
    ////////////////////////////////////////////////////////////////////////////////
    Real_t width;
    if(JSONHelpers::readValue(jParams, width, "Width")) {
        setWidth(width);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
Real_t CapsuleObject<N, Real_t>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const {
    auto ppos = this->invTransform(ppos0);
    VecN pa   = ppos - m_Start;
    VecN ba   = m_End - m_Start;

    Real_t h = MathHelpers::clamp(glm::dot(pa, ba) / glm::dot(ba, ba), Real_t(0.0), Real_t(1.0));
    Real_t d = this->m_UniformScale * (glm::length(pa - ba * h) - m_Radius);
    return bNegativeInside ? d : -d;
}

////////////////////////////////////////////////////////////////////////////////
template<Int N, class Real_t>
void CapsuleObject<N, Real_t>::parseParameters(const JParams& jParams) {
    GeometryObject<N, Real_t>::parseParameters(jParams);
    ////////////////////////////////////////////////////////////////////////////////
    if(Real_t radius; JSONHelpers::readValue(jParams, radius, "Radius")) {
        setRadius(radius);
    }
    ////////////////////////////////////////////////////////////////////////////////
    m_Start[0] = Real_t(-1);
    m_End[0]   = Real_t(1);
    JSONHelpers::readValue(jParams, m_Start[0], "Start");
    JSONHelpers::readValue(jParams, m_End[0],   "End");
    m_Start[0] += m_Radius;
    m_End[0]   -= m_Radius;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
Real_t EllipsoidObject<N, Real_t>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const {
    auto ppos = this->invTransform(ppos0);
    for(Int d = 0; d < N; ++d) {
        ppos[d] /= m_RadiusRatio[d];
    }
    Real_t d = this->m_UniformScale * glm::compMin(m_RadiusRatio) * (glm::length(ppos) - Real_t(1.0));
    return bNegativeInside ? d : -d;
}

////////////////////////////////////////////////////////////////////////////////
template<Int N, class Real_t>
void EllipsoidObject<N, Real_t>::parseParameters(const JParams& jParams) {
    GeometryObject<N, Real_t>::parseParameters(jParams);
    ////////////////////////////////////////////////////////////////////////////////
    VecN radiusRatio;
    if(JSONHelpers::readVector(jParams, radiusRatio, "RadiusRatio")) {
        setRadiusRatio(radiusRatio);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// sign distance field for triangle mesh
template<class Real_t>
void computeSDFMesh(const StdVT<Vec3ui>& faces, const StdVT_Vec3<Real_t>& vertices, const Vec3<Real_t>& origin, Real_t cellSize,
                    UInt ni, UInt nj, UInt nk, Array<3, Real_t>& SDF, Int exactBand = 1) {
    NT_REQUIRE(ni > 0 && nj > 0 && nk > 0);

    SDF.resize(ni, nj, nk);
    SDF.assign(Real_t(ni + nj + nk) * cellSize); // upper bound on distance
    Array3ui closest_tri(ni, nj, nk, 0xffffffff);

    // intersection_count(i,j,k) is # of tri intersections in (i-1,i]x{j}x{k}
    // we begin by initializing distances near the mesh, and figuring out intersection counts
    Array3ui intersectionCount(ni, nj, nk, 0u);

    for(UInt face = 0, faceEnd = static_cast<UInt>(faces.size()); face < faceEnd; ++face) {
        UInt p = faces[face][0];
        UInt q = faces[face][1];
        UInt r = faces[face][2];

        // coordinates in grid to high precision
        Vec3<Real_t> fp = (vertices[p] - origin) / cellSize;
        Vec3<Real_t> fq = (vertices[q] - origin) / cellSize;
        Vec3<Real_t> fr = (vertices[r] - origin) / cellSize;

        // do distances nearby
        Int i0 = MathHelpers::clamp(static_cast<Int>(MathHelpers::min(fp[0], fq[0], fr[0])) - exactBand, 0, static_cast<Int>(ni - 1));
        Int i1 = MathHelpers::clamp(static_cast<Int>(MathHelpers::max(fp[0], fq[0], fr[0])) + exactBand + 1, 0, static_cast<Int>(ni - 1));
        Int j0 = MathHelpers::clamp(static_cast<Int>(MathHelpers::min(fp[1], fq[1], fr[1])) - exactBand, 0, static_cast<Int>(nj - 1));
        Int j1 = MathHelpers::clamp(static_cast<Int>(MathHelpers::max(fp[1], fq[1], fr[1])) + exactBand + 1, 0, static_cast<Int>(nj - 1));
        Int k0 = MathHelpers::clamp(static_cast<Int>(MathHelpers::min(fp[2], fq[2], fr[2])) - exactBand, 0, static_cast<Int>(nk - 1));
        Int k1 = MathHelpers::clamp(static_cast<Int>(MathHelpers::max(fp[2], fq[2], fr[2])) + exactBand + 1, 0, static_cast<Int>(nk - 1));

        ParallelExec::run<Int>(i0, i1 + 1, j0, j1 + 1, k0, k1 + 1,
                               [&](Int i, Int j, Int k) {
                                   Vec3<Real_t> gx = Vec3<Real_t>(i, j, k) * cellSize + origin;
                                   Real_t d        = GeometryHelpers::point_triangle_distance(gx, vertices[p], vertices[q], vertices[r]);

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
                Real_t a, b, c;

                if(GeometryHelpers::point_in_triangle_2d(static_cast<Real_t>(j), static_cast<Real_t>(k), fp[1], fp[2], fq[1], fq[2], fr[1], fr[2], a, b, c)) {
                    // intersection i coordinate
                    Real_t fi = a * fp[0] + b * fq[0] + c * fr[0];

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
    } // end loop face

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
    ParallelExec::run<UInt>(0, nk,
                            [&](UInt k) {
                                for(UInt j = 0; j < nj; ++j) {
                                    UInt total_count = 0;

                                    for(UInt i = 0; i < ni; ++i) {
                                        total_count += intersectionCount(i, j, k);

                                        if(total_count & 1) {             // if parity of intersections so far is odd,
                                            SDF(i, j, k) = -SDF(i, j, k); // we are inside the mesh
                                        }
                                    }
                                }
                            });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
Real_t TriMeshObject<N, Real_t>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const {
    if constexpr(N == 2) {
        NT_UNUSED(ppos0);
        NT_UNUSED(bNegativeInside);
        return 0;
    } else {
        NT_REQUIRE(m_bSDFGenerated);
        auto   ppos    = this->invTransform(ppos0);
        auto   gridPos = m_Grid3D.getGridCoordinate(ppos);
        Real_t d       = this->m_UniformScale * ArrayHelpers::interpolateValueLinear(gridPos, m_SDFData);
        return bNegativeInside ? d : -d;
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void TriMeshObject<N, Real_t>::computeSDF() {
    if constexpr(N == 3) {
        ////////////////////////////////////////////////////////////////////////////////
        // Load mesh
        MeshLoader meshLoader;
        NT_REQUIRE(meshLoader.loadMesh(m_TriMeshFile));
        meshLoader.scaleToBox();

        ////////////////////////////////////////////////////////////////////////////////
        // Init grid and array of data
        m_Grid3D.setGrid(VecN(meshLoader.getAABBMin()) - VecN(Real_t(3.0) * m_Step),
                         VecN(meshLoader.getAABBMax()) + VecN(Real_t(3.0) * m_Step),
                         m_Step);

        StdVT_Vec3<Real_t> vertexList(meshLoader.getNVertices());
        StdVT<Vec3ui>      faceList(meshLoader.getNFaces());

        std::memcpy(vertexList.data(), meshLoader.getVertices().data(), meshLoader.getVertices().size() * sizeof(Real_t));
        std::memcpy(faceList.data(),   meshLoader.getFaces().data(),    meshLoader.getFaces().size() * sizeof(UInt));

        ////////////////////////////////////////////////////////////////////////////////
        // Compute SDF data
        computeSDFMesh(faceList, vertexList, VecN(meshLoader.getAABBMin()),
                       m_Step, m_Grid3D.getNCells()[0], m_Grid3D.getNCells()[1], m_Grid3D.getNCells()[2], m_SDFData);
        m_bSDFGenerated = true;
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void TriMeshObject<N, Real_t>::parseParameters(const JParams& jParams) {
    GeometryObject<N, Real_t>::parseParameters(jParams);
    ////////////////////////////////////////////////////////////////////////////////
    NT_REQUIRE(JSONHelpers::readValue(jParams, meshFile(), "MeshFile"));
    JSONHelpers::readValue(jParams, sdfStep(), "SDFStep");
    computeSDF();
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
Real_t CSGObject<N, Real_t>::signedDistance(const VecN& ppos0, bool bNegativeInside /*= true*/) const {
    if(m_Objects.size() == 0) {
        return HugeReal();
    }

    auto   ppos = domainDeform(this->transform(ppos0));
    Real_t sd   = m_Objects[0].obj->signedDistance(ppos, bNegativeInside);

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
template<Int N, class Real_t>
VecX<N, Real_t> CSGObject<N, Real_t>::domainDeform(const VecN& ppos) const {
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
template<Int N, class Real_t>
VecX<N, Real_t> CSGObject<N, Real_t>::twist(const VecN& ppos) const {
    NT_UNUSED(ppos);
    //Real_t         c = cos(Real_t(5.0) * ppos.z);
    //Real_t         s = sin(Real_t(5.0) * ppos.z);
    //Mat2x2<Real_t> m = Mat2x2<Real_t>(c, -s, s, c);
    //Mat2x2<Real_t> m = Mat2x2<Real_t>(c, s, -s, c);

    //return VecN(m * VecN(ppos[0], ppos[1]), ppos.z);
    return VecN(0);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
VecX<N, Real_t> CSGObject<N, Real_t>::cheapBend(const VecN& ppos) const {
    Real_t c = cos(Real_t(0.5) * ppos.y);
    Real_t s = sin(Real_t(0.5) * ppos.y);
    //Mat2x2<Real_t> m = Mat2x2<Real_t>(c, -s, s, c);
    Mat2x2<Real_t> m = Mat2x2<Real_t>(c, s, -s, c);

    //return VecN(m * VecN(ppos[0], ppos[1]), ppos.z);
    return VecN(0);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(GeometryObject)
NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(BoxObject)
NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(SphereObject)
NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(TorusObject)
NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(Torus28Object)
NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(Torus2InfObject)
NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(Torus88Object)
NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(TorusInfInfObject)
NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(CylinderObject)
NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(ConeObject)
NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(PlaneObject)
NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(TriangleObject)
NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(HexagonObject)
NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(TriangularPrismObject)
NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(HexagonalPrismObject)
NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(CapsuleObject)
NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(EllipsoidObject)
NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(TriMeshObject)
NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(CSGObject)

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace NTCodeBase
