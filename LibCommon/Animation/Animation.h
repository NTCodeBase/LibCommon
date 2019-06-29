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
#include <LibCommon/Animation/CubicSpline.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// utility to hack glm error
namespace glm {
template<class T> NTCodeBase::Mat3x3<T> rotate(const NTCodeBase::Mat3x3<T>& m, T angle, const NTCodeBase::Vec2<T>&) {
    return glm::rotate(m, angle);
}
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NTCodeBase {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
struct KeyFrame {
    ////////////////////////////////////////////////////////////////////////////////
    NT_TYPE_ALIAS
    ////////////////////////////////////////////////////////////////////////////////
    KeyFrame() = default;
    KeyFrame(UInt frame_, const VecN& translation_) : frame(frame_), translation(translation_) {}
    KeyFrame(UInt frame_, const VecNp1& rotation_) : frame(frame_), rotation(rotation_) {}
    KeyFrame(UInt frame_, Real_t scale_) : frame(frame_), uniformScale(scale_), invScale(Real_t(1.0) / scale_) {}
    KeyFrame(UInt frame_, const VecN& translation_, const VecNp1& rotation_, Real_t scale_ = Real_t(1.0))
        : frame(frame_), translation(translation_), rotation(rotation_), uniformScale(scale_), invScale(Real_t(1.0) / scale_) {}
    ////////////////////////////////////////////////////////////////////////////////
    UInt   frame        = 0;
    VecN   translation  = VecN(0);
    VecNp1 rotation     = VecNp1(VecN(1), 0);
    Real_t uniformScale = Real_t(1.0);
    Real_t invScale     = Real_t(1.0);
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class RigidBodyAnimation {
    ////////////////////////////////////////////////////////////////////////////////
    NT_TYPE_ALIAS
    ////////////////////////////////////////////////////////////////////////////////
public:
    RigidBodyAnimation() = default;
    ////////////////////////////////////////////////////////////////////////////////
    auto isActive(UInt frame) { return (m_KeyFrames.size() > 0) && (frame >= m_StartFrame) && (m_bPeriodic || frame <= m_EndFrame); }
    auto doneAnimation(UInt frame) const { return (frame > m_EndFrame && !m_bPeriodic); }
    auto& periodic() { return m_bPeriodic; }
    auto& keyFrames() { return m_KeyFrames; }
    auto nKeyFrames() const { return static_cast<UInt>(m_KeyFrames.size()); }
    ////////////////////////////////////////////////////////////////////////////////
    auto& addKeyFrame() { return m_KeyFrames.emplace_back(KeyFrame<N, Real_t>()); }
    void addKeyFrame(const KeyFrame<N, Real_t>& keyFrame) { m_KeyFrames.push_back(keyFrame); }
    void addKeyFrame(UInt frame, const VecN& translation) { m_KeyFrames.emplace_back(KeyFrame<N, Real_t>(frame, translation)); }
    void addKeyFrame(UInt frame, const VecNp1& rotation) { m_KeyFrames.emplace_back(KeyFrame<N, Real_t>(frame, rotation)); }
    void addKeyFrame(UInt frame, const VecN& translation, const VecNp1& rotation) { m_KeyFrames.emplace_back(KeyFrame<N, Real_t>(frame, translation, rotation)); }
    ////////////////////////////////////////////////////////////////////////////////
    void       makeReady(bool bCubicIntTranslation = true, bool bCubicIntRotation = true);
    MatNp1xNp1 getInvTransformation(UInt frame, Real_t frameFraction = Real_t(0));
    MatNp1xNp1 getInactiveTransformationMatrix(UInt frame);
    ////////////////////////////////////////////////////////////////////////////////
    virtual MatNp1xNp1 getTransformationMatrix(UInt frame, Real_t frameFraction = Real_t(0));

protected:
    StdVT<KeyFrame<N, Real_t>> m_KeyFrames;
    CubicSpline<Real_t>        m_TranslationInterpolator[N];
    CubicSpline<Real_t>        m_RotationInterpolator[N + 1];

    MatNp1xNp1 m_StartFrameTransformationMatrix = MatNp1xNp1(1);
    MatNp1xNp1 m_EndFrameTransformationMatrix   = MatNp1xNp1(1);

    UInt m_StartFrame = 0;
    UInt m_EndFrame   = 1u;
    UInt m_FrameRange = 1u;
    bool m_bReady     = false;
    bool m_bPeriodic  = false;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class Animation : public RigidBodyAnimation<N, Real_t> {
    ////////////////////////////////////////////////////////////////////////////////
    NT_TYPE_ALIAS
    ////////////////////////////////////////////////////////////////////////////////
public:
    Animation() = default;
    ////////////////////////////////////////////////////////////////////////////////
    void addKeyFrame(UInt frame, Real_t scale) { this->m_KeyFrames.emplace_back(KeyFrame<N, Real_t>(frame, scale)); }
    void addKeyFrame(UInt frame, const VecN& translation, const VecNp1& rotation, Real_t scale) {
        this->m_KeyFrames.emplace_back(KeyFrame<N, Real_t>(frame, translation, rotation, scale));
    }

    ////////////////////////////////////////////////////////////////////////////////
    void   makeReady(bool bCubicIntTranslation = true, bool bCubicIntRotation = true, bool bCubicIntScale = true);
    Real_t getUniformScale(UInt frame, Real_t frameFraction = Real_t(0));
    ////////////////////////////////////////////////////////////////////////////////
    virtual MatNp1xNp1 getTransformationMatrix(UInt frame, Real_t frameFraction = Real_t(0)) override;

protected:
    CubicSpline<Real_t> m_ScaleInterpolator;
};
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace NTCodeBase
