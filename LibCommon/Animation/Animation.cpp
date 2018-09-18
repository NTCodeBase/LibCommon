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
void RigidBodyAnimation<N, RealType>::setAnimationRange(UInt startFrame, UInt endFrame)
{
    m_StartFrame = startFrame;
    if(endFrame > 0) {
        m_EndFrame   = endFrame;
        m_FrameRange = m_EndFrame - m_StartFrame;
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void RigidBodyAnimation<N, RealType>::makeReady(bool bCubicIntTranslation, bool bCubicIntRotation)
{
    if(m_KeyFrames.size() == 0) {
        return;
    }
    ////////////////////////////////////////////////////////////////////////////////
    // add more 1 point if there are only 2 key frames
    if(m_KeyFrames.size() == 2) {
        auto& keyFrame  = m_KeyFrames.emplace_back(KeyFrame<N, RealType>());
        auto& keyFrame0 = m_KeyFrames[0];
        auto& keyFrame1 = m_KeyFrames[1];
        keyFrame.frame        = (keyFrame0.frame + keyFrame1.frame) / 2;
        keyFrame.translation  = (keyFrame0.translation + keyFrame1.translation) * RealType(0.5);
        keyFrame.rotation     = (keyFrame0.rotation + keyFrame1.rotation) * RealType(0.5);
        keyFrame.uniformScale = (keyFrame0.uniformScale + keyFrame1.uniformScale) * RealType(0.5);
        keyFrame.invScale     = RealType(1) / keyFrame.uniformScale;
        std::swap(keyFrame, keyFrame1);
    }

    ////////////////////////////////////////////////////////////////////////////////
    StdVT<RealType> frames;
    StdVT<RealType> translations[N];
    StdVT<RealType> rotations[N + 1];

    frames.reserve(nKeyFrames());
    for(Int d = 0; d < N; ++d) {
        translations[d].reserve(nKeyFrames());
        rotations[d].reserve(nKeyFrames());
    }
    rotations[N].reserve(nKeyFrames());
    ////////////////////////////////////////////////////////////////////////////////
    UInt maxFrame = 0;
    for(const auto& keyFrame : m_KeyFrames) {
        frames.push_back(static_cast<RealType>(keyFrame.frame));
        for(Int d = 0; d < N; ++d) {
            translations[d].push_back(keyFrame.translation[d]);
            rotations[d].push_back(keyFrame.rotation[d]);
        }
        rotations[N].push_back(keyFrame.rotation[N]);
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
    ////////////////////////////////////////////////////////////////////////////////
    __NT_REQUIRE(m_EndFrame > m_StartFrame);
    m_FrameRange = m_EndFrame - m_StartFrame;
    m_bReady     = true;
    ////////////////////////////////////////////////////////////////////////////////
    m_StartFrameTransformationMatrix = getTransformationMatrix(m_StartFrame);
    m_EndFrameTransformationMatrix   = getTransformationMatrix(m_EndFrame);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
MatXxX<N + 1, RealType> RigidBodyAnimation<N, RealType>::getInactiveTransformationMatrix(UInt frame)
{
    if(frame < m_StartFrame) {
        return m_StartFrameTransformationMatrix;
    } else if(frame > m_EndFrame) {
        return m_EndFrameTransformationMatrix;
    } else {
        return getTransformationMatrix(frame);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
MatXxX<N + 1, RealType> RigidBodyAnimation<N, RealType>::getTransformationMatrix(UInt frame, RealType frameFraction /*= RealType(0)*/)
{
    if(nKeyFrames() == 0) {
        return MatNp1xNp1(1);
    }
    ////////////////////////////////////////////////////////////////////////////////
    __NT_REQUIRE(m_bReady)
    VecN translation;
    VecNp1 rotation;

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

    MatNp1xNp1 transMatrix(1);
    if(rotation[N] > 0) {
        transMatrix = glm::rotate(transMatrix, rotation[N], VecN(rotation));
    }
    transMatrix = glm::translate(transMatrix, translation);
    return transMatrix;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
MatXxX<N + 1, RealType> RigidBodyAnimation<N, RealType>::getInvTransformation(UInt frame, RealType frameFraction /*= RealType(0)*/)
{
    if(nKeyFrames() == 0) {
        return MatNp1xNp1(1);
    }
    ////////////////////////////////////////////////////////////////////////////////
    __NT_REQUIRE(m_bReady);
    return glm::inverse(getTransformationMatrix(frame, frameFraction));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
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
template<Int N, class RealType>
MatXxX<N + 1, RealType> Animation<N, RealType>::getTransformationMatrix(UInt frame, RealType frameFraction /*= RealType(0)*/)
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
    if(rotation[N] > 0) {
        transMatrix = glm::rotate(transMatrix, rotation[N], VecN(rotation));
    }
    transMatrix = glm::translate(transMatrix, translation);
    return transMatrix;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(RigidBodyAnimation)
