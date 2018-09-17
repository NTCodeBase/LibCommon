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
    KeyFrame(UInt frame_, const VecX<N + 1, RealType>& rotation_) : frame(frame_), rotation(rotation_) {}
    KeyFrame(UInt frame_, RealType scale_) : frame(frame_), uniformScale(scale_), invScale(RealType(1.0) / scale_) {}
    KeyFrame(UInt frame_, const VecN& translation_, const VecX<N + 1, RealType>& rotation_, RealType scale_)
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
class Animation
{
    ////////////////////////////////////////////////////////////////////////////////
    __NT_TYPE_ALIASING
    ////////////////////////////////////////////////////////////////////////////////
public:
    Animation() { m_KeyFrames.emplace_back(KeyFrame<N, RealType>()); }
    ////////////////////////////////////////////////////////////////////////////////
    void addKeyFrame(const KeyFrame<N, RealType>& keyFrame);
    void addKeyFrame(UInt frame, const VecN& translation);
    void addKeyFrame(UInt frame, const VecX<N + 1, RealType>& rotation);
    void addKeyFrame(UInt frame, RealType scale);
    void addKeyFrame(UInt frame, const VecN& translation, const VecX<N + 1, RealType>& rotation, RealType scale);

    UInt  nKeyFrames() const { return static_cast<UInt>(m_KeyFrames.size()); }
    auto& keyFrames() { return m_KeyFrames; }
    void  setPeriodic(bool bPeriodic, UInt startFrame = 0) { m_bPeriodic = bPeriodic; m_StartFrame = startFrame; }

    void                    makeReady(bool bCubicIntTranslation = true, bool bCubicIntRotation = true, bool bCubicIntScale = true);
    void                    getTransformation(VecN& translation, VecX<N + 1, RealType>& rotation, RealType& scale, UInt frame, RealType fraction = RealType(0));
    MatXxX<N + 1, RealType> getTransformation(UInt frame, RealType fraction    = RealType(0));
    MatXxX<N + 1, RealType> getInvTransformation(UInt frame, RealType fraction = RealType(0));
    RealType                getUniformScale(UInt frame, RealType fraction      = RealType(0));

private:
    StdVT<KeyFrame<N, RealType>> m_KeyFrames;
    CubicSpline<RealType>        m_TranslationSpline[N];
    CubicSpline<RealType>        m_RotationSpline[N + 1];
    CubicSpline<RealType>        m_ScaleSpline;

    UInt m_StartFrame = 0;
    UInt m_MaxFrame   = 0;
    bool m_bReady     = false;
    bool m_bPeriodic  = false;
};
