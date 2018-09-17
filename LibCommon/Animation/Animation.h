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
#include <LibCommon/Animation/CubicSpline.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// utility to hack glm error
namespace glm
{
template<class T> Mat3x3<T> rotate(const Mat3x3<T>& m, T angle, const Vec2<T>&) { return glm::rotate(m, angle); }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
struct KeyFrame
{
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
    KeyFrame() = default;
    KeyFrame(UInt frame_, const VecN& translation_) : frame(frame_), translation(translation_) {}
    KeyFrame(UInt frame_, const VecNp1& rotation_) : frame(frame_), rotation(rotation_) {}
    KeyFrame(UInt frame_, RealType scale_) : frame(frame_), uniformScale(scale_), invScale(RealType(1.0) / scale_) {}
    KeyFrame(UInt frame_, const VecN& translation_, const VecNp1& rotation_, RealType scale_ = RealType(1.0))
        : frame(frame_), translation(translation_), rotation(rotation_), uniformScale(scale_), invScale(RealType(1.0) / scale_) {}
    ////////////////////////////////////////////////////////////////////////////////
    UInt     frame        = 0;
    VecN     translation  = VecN(0);
    VecNp1   rotation     = VecNp1(VecN(1), 0);
    RealType uniformScale = RealType(1.0);
    RealType invScale     = RealType(1.0);
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
class RigidBodyAnimation
{
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
public:
    RigidBodyAnimation() = default;
    ////////////////////////////////////////////////////////////////////////////////
    auto isActive(UInt frame) { return (m_KeyFrames.size() > 0) && (frame >= m_StartFrame) && (m_bPeriodic || frame <= m_EndFrame); }
    auto nKeyFrames() const { return static_cast<UInt>(m_KeyFrames.size()); }
    auto& keyFrames() { return m_KeyFrames; }
    void setPeriodic(bool bPeriodic) { m_bPeriodic = bPeriodic; }
    void setAnimationRange(UInt startFrame, UInt endFrame = 0);

    ////////////////////////////////////////////////////////////////////////////////
    auto& addKeyFrame() { return m_KeyFrames.emplace_back(KeyFrame<N, RealType>()); }
    void addKeyFrame(const KeyFrame<N, RealType>& keyFrame) { m_KeyFrames.push_back(keyFrame); }
    void addKeyFrame(UInt frame, const VecN& translation) { m_KeyFrames.emplace_back(KeyFrame<N, RealType>(frame, translation)); }
    void addKeyFrame(UInt frame, const VecNp1& rotation) { m_KeyFrames.emplace_back(KeyFrame<N, RealType>(frame, rotation)); }
    void addKeyFrame(UInt frame, const VecN& translation, const VecNp1& rotation) { m_KeyFrames.emplace_back(KeyFrame<N, RealType>(frame, translation, rotation)); }
    ////////////////////////////////////////////////////////////////////////////////
    void       makeReady(bool bCubicIntTranslation = true, bool bCubicIntRotation = true);
    MatNp1xNp1 getInvTransformation(UInt frame, RealType frameFraction = RealType(0));
    ////////////////////////////////////////////////////////////////////////////////
    virtual MatNp1xNp1 getTransformationMatrix(UInt frame, RealType frameFraction = RealType(0));

protected:
    StdVT<KeyFrame<N, RealType>> m_KeyFrames;
    CubicSpline<RealType>        m_TranslationInterpolator[N];
    CubicSpline<RealType>        m_RotationInterpolator[N + 1];

    UInt m_StartFrame = 0;
    UInt m_EndFrame   = 1u;
    UInt m_FrameRange = 1u;
    bool m_bReady     = false;
    bool m_bPeriodic  = false;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
class Animation : public RigidBodyAnimation<N, RealType>
{
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
public:
    Animation() = default;
    ////////////////////////////////////////////////////////////////////////////////
    void addKeyFrame(UInt frame, RealType scale) { m_KeyFrames.emplace_back(KeyFrame<N, RealType>(frame, scale)); }
    void addKeyFrame(UInt frame, const VecN& translation, const VecNp1& rotation, RealType scale)
    {
        m_KeyFrames.emplace_back(KeyFrame<N, RealType>(frame, translation, rotation, scale));
    }

    ////////////////////////////////////////////////////////////////////////////////
    void     makeReady(bool bCubicIntTranslation = true, bool bCubicIntRotation = true, bool bCubicIntScale = true);
    RealType getUniformScale(UInt frame, RealType frameFraction = RealType(0));
    ////////////////////////////////////////////////////////////////////////////////
    virtual MatNp1xNp1 getTransformationMatrix(UInt frame, RealType frameFraction = RealType(0)) override;

protected:
    CubicSpline<RealType> m_ScaleInterpolator;
};
