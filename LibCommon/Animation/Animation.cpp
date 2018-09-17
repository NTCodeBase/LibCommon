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

#include <LibCommon/Animation/Animation.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Animation<N, RealType>::setAnimationRange(UInt startFrame, UInt endFrame)
{
    m_StartFrame = startFrame;
    if(endFrame > 0) {
        m_EndFrame   = endFrame;
        m_FrameRange = m_EndFrame - m_StartFrame;
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Animation<N, RealType>::makeReady(bool bCubicIntTranslation, bool bCubicIntRotation, bool bCubicIntScale)
{
    if(m_KeyFrames.size() == 0) {
        return;
    }
    ////////////////////////////////////////////////////////////////////////////////
    StdVT<RealType> frames;
    StdVT<RealType> translations[N];
    StdVT<RealType> rotations[N + 1];
    StdVT<RealType> scales;

    frames.reserve(nKeyFrames());
    for(Int d = 0; d < N; ++d) {
        translations[d].reserve(nKeyFrames());
        rotations[d].reserve(nKeyFrames());
    }
    rotations[N].reserve(nKeyFrames());
    scales.reserve(nKeyFrames());
    ////////////////////////////////////////////////////////////////////////////////
    UInt maxFrame = 0;
    for(const auto& keyFrame : m_KeyFrames) {
        frames.push_back(static_cast<RealType>(keyFrame.frame));
        for(Int d = 0; d < N; ++d) {
            translations[d].push_back(keyFrame.translation[d]);
            rotations[d].push_back(keyFrame.rotation[d]);
        }
        rotations[N].push_back(keyFrame.rotation[N]);
        scales.push_back(keyFrame.uniformScale);
        ////////////////////////////////////////////////////////////////////////////////
        if(maxFrame < keyFrame.frame) {
            maxFrame = keyFrame.frame;
        }
    }
    if(m_EndFrame == 0) {
        // if end frame has not been set, set it to the latest key frame
        m_EndFrame = maxFrame;
    }
    ////////////////////////////////////////////////////////////////////////////////
    for(Int d = 0; d < N; ++d) {
        m_TranslationInterpolator[d].setBoundary(CubicSpline<RealType>::BDType::FirstOrder, 0, CubicSpline<RealType>::BDType::FirstOrder, 0);
        m_TranslationInterpolator[d].setPoints(frames, translations[d], bCubicIntTranslation);

        m_RotationInterpolator[d].setBoundary(CubicSpline<RealType>::BDType::FirstOrder, 0, CubicSpline<RealType>::BDType::FirstOrder, 0);
        m_RotationInterpolator[d].setPoints(frames, rotations[d], bCubicIntRotation);
    }
    m_RotationInterpolator[N].setBoundary(CubicSpline<RealType>::BDType::FirstOrder, 0, CubicSpline<RealType>::BDType::FirstOrder, 0);
    m_RotationInterpolator[N].setPoints(frames, rotations[N], bCubicIntRotation);
    m_ScaleInterpolator.setBoundary(CubicSpline<RealType>::BDType::FirstOrder, 0, CubicSpline<RealType>::BDType::FirstOrder, 0);
    m_ScaleInterpolator.setPoints(frames, scales, bCubicIntScale);
    ////////////////////////////////////////////////////////////////////////////////
    __NT_REQUIRE(m_EndFrame > m_StartFrame);
    m_FrameRange = m_EndFrame - m_StartFrame;
    m_bReady     = true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Animation<N, RealType>::getTransformation(VecN& translation, VecNp1& rotation, RealType& scale, UInt frame, RealType frameFraction /*= RealType(0)*/)
{
    if(nKeyFrames() == 0) {
        translation = VecN(0);
        rotation    = VecNp1(0);
        scale       = RealType(1);
        return;
    }
    ////////////////////////////////////////////////////////////////////////////////
    __NT_REQUIRE(m_bReady)
    if(frame < m_StartFrame) {
        frame = m_StartFrame;
    } else if(frame >= m_EndFrame) {
        if(m_bPeriodic) {
            frame = ((frame - m_StartFrame) % m_FrameRange) + m_StartFrame;
        } else {
            frame         = m_EndFrame;
            frameFraction = RealType(0);
        }
    }

    RealType x = static_cast<RealType>(frame) + frameFraction;
    for(Int d = 0; d < N; ++d) {
        translation[d] = m_TranslationInterpolator[d](x);
        rotation[d]    = m_RotationInterpolator[d](x);
    }
    rotation[N] = m_RotationInterpolator[N](x);
    scale       = m_ScaleInterpolator(x);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
MatXxX<N + 1, RealType> Animation<N, RealType>::getTransformation(UInt frame, RealType frameFraction /*= RealType(0)*/)
{
    if(nKeyFrames() == 0) {
        return MatNp1xNp1(1);
    }
    ////////////////////////////////////////////////////////////////////////////////
    __NT_REQUIRE(m_bReady)
    VecN translation;
    VecNp1   rotation;
    RealType scale;

    if(frame < m_StartFrame) {
        frame = m_StartFrame;
    } else if(frame >= m_EndFrame) {
        if(m_bPeriodic) {
            frame = ((frame - m_StartFrame) % m_FrameRange) + m_StartFrame;
        } else {
            frame         = m_EndFrame;
            frameFraction = RealType(0);
        }
    }

    RealType x = static_cast<RealType>(frame) + frameFraction;
    for(Int d = 0; d < N; ++d) {
        translation[d] = m_TranslationInterpolator[d](x);
        rotation[d]    = m_RotationInterpolator[d](x);
    }
    rotation[N] = m_RotationInterpolator[N](x);
    scale       = m_ScaleInterpolator(x);

    MatNp1xNp1 transMatrix(1);
    transMatrix = glm::scale(transMatrix, VecN(scale));
    transMatrix = glm::rotate(transMatrix, rotation[N], VecN(rotation));
    transMatrix = glm::translate(transMatrix, translation);
    return transMatrix;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
MatXxX<N + 1, RealType> Animation<N, RealType>::getInvTransformation(UInt frame, RealType frameFraction /*= RealType(0)*/)
{
    if(nKeyFrames() == 0) {
        return MatNp1xNp1(1);
    }
    ////////////////////////////////////////////////////////////////////////////////
    __NT_REQUIRE(m_bReady);
    return glm::inverse(getTransformation(frame, frameFraction));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
RealType Animation<N, RealType>::getUniformScale(UInt frame, RealType frameFraction /*= RealType(0)*/)
{
    if(nKeyFrames() == 0) {
        return RealType(1);
    }
    ////////////////////////////////////////////////////////////////////////////////
    if(frame < m_StartFrame) {
        frame = m_StartFrame;
    } else if(frame >= m_EndFrame) {
        if(m_bPeriodic) {
            frame = ((frame - m_StartFrame) % m_FrameRange) + m_StartFrame;
        } else {
            frame         = m_EndFrame;
            frameFraction = RealType(0);
        }
    }
    RealType x = static_cast<RealType>(frame) + frameFraction;
    return m_ScaleInterpolator(x);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(Animation)
