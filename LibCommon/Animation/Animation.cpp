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
void Animation<N, RealType>::addKeyFrame(const KeyFrame<N, RealType>& keyFrame)
{
    if(keyFrame.frame == 0) {
        m_KeyFrames[0] = keyFrame;
    } else {
        m_KeyFrames.push_back(keyFrame);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Animation<N, RealType>::addKeyFrame(UInt frame, const VecN& translation)
{
    if(frame == 0) {
        m_KeyFrames[0].translation = translation;
    } else {
        m_KeyFrames.emplace_back(KeyFrame<N, RealType>(frame, translation));
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Animation<N, RealType>::addKeyFrame(UInt frame, const VecNp1& rotation)
{
    if(frame == 0) {
        m_KeyFrames[0].rotation = rotation;
    } else {
        m_KeyFrames.emplace_back(KeyFrame<N, RealType>(frame, rotation));
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Animation<N, RealType>::addKeyFrame(UInt frame, RealType scale)
{
    if(frame == 0) {
        m_KeyFrames[0].uniformScale = scale;
    } else {
        m_KeyFrames.emplace_back(KeyFrame<N, RealType>(frame, scale));
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Animation<N, RealType>::addKeyFrame(UInt frame, const VecN& translation, const VecNp1& rotation, RealType scale)
{
    if(frame == 0) {
        m_KeyFrames[0].translation  = translation;
        m_KeyFrames[0].rotation     = rotation;
        m_KeyFrames[0].uniformScale = scale;
    } else {
        m_KeyFrames.emplace_back(KeyFrame<N, RealType>(frame, translation, rotation, scale));
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Animation<N, RealType>::makeReady(bool bCubicIntTranslation, bool bCubicIntRotation, bool bCubicIntScale)
{
    StdVT<RealType> frames;
    StdVT<RealType> translations[N];
    StdVT<RealType> rotations[N + 1];
    StdVT<RealType> scales;

    for(Int d = 0; d < N; ++d) {
        translations[d].reserve(nKeyFrames());
        rotations[d].reserve(nKeyFrames());
    }
    rotations[N].reserve(nKeyFrames());
    frames.reserve(nKeyFrames());
    scales.reserve(nKeyFrames());
    ////////////////////////////////////////////////////////////////////////////////
    for(const auto& keyFrame : m_KeyFrames) {
        if(m_bPeriodic && keyFrame.frame < m_StartFrame) { // ignore frame before start frame
            continue;
        }
        m_EndFrame = (m_EndFrame < keyFrame.frame) ? keyFrame.frame : m_EndFrame; // set end frame to the latest frame if applicable

        for(Int d = 0; d < N; ++d) {
            translations[d].push_back(keyFrame.translation[d]);
            rotations[d].push_back(keyFrame.rotation[d]);
        }
        rotations[N].push_back(keyFrame.rotation[N]);
        scales.push_back(keyFrame.uniformScale);
        frames.push_back(static_cast<RealType>(keyFrame.frame));
    }

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
    m_bReady = true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Animation<N, RealType>::getTransformation(VecN& translation, VecNp1& rotation, RealType& scale, UInt frame, RealType frameFraction /*= RealType(0)*/)
{
    if(m_KeyFrames.size() == 1 || (m_bPeriodic && frame < m_StartFrame)) {
        translation = m_KeyFrames[0].translation;
        rotation    = m_KeyFrames[0].rotation;
        scale       = m_KeyFrames[0].uniformScale;
        return;
    }
    ////////////////////////////////////////////////////////////////////////////////
    __NT_REQUIRE(m_bReady)
    if(m_bPeriodic && frame >= m_EndFrame) {
        frame = ((frame - m_StartFrame) % (m_EndFrame - m_StartFrame)) + m_StartFrame;
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
    if(m_KeyFrames.size() == 1 || (m_bPeriodic && frame < m_StartFrame)) {
        MatNp1xNp1 translationMatrix = glm::translate(MatNp1xNp1(1.0), m_KeyFrames[0].translation);
        MatNp1xNp1 rotationMatrix    = glm::rotate(MatNp1xNp1(1.0), m_KeyFrames[0].rotation[N], VecN(m_KeyFrames[0].rotation));
        MatNp1xNp1 scaleMatrix       = glm::scale(MatNp1xNp1(1.0), VecN(m_KeyFrames[0].uniformScale));
        return (translationMatrix * rotationMatrix * scaleMatrix);
    }
    ////////////////////////////////////////////////////////////////////////////////
    __NT_REQUIRE(m_bReady)
    VecN translation;
    VecNp1   rotation;
    RealType scale;

    if(m_bPeriodic && frame >= m_EndFrame) {
        frame = ((frame - m_StartFrame) % (m_EndFrame - m_StartFrame)) + m_StartFrame;
    }
    RealType x = static_cast<RealType>(frame) + frameFraction;

    for(Int d = 0; d < N; ++d) {
        translation[d] = m_TranslationInterpolator[d](x);
        rotation[d]    = m_RotationInterpolator[d](x);
    }
    rotation[N] = m_RotationInterpolator[N](x);
    scale       = m_ScaleInterpolator(x);

    MatNp1xNp1 translationMatrix = glm::translate(MatNp1xNp1(1.0), translation);
    MatNp1xNp1 rotationMatrix    = glm::rotate(MatNp1xNp1(1.0), rotation[N], VecN(rotation));
    MatNp1xNp1 scaleMatrix       = glm::scale(MatNp1xNp1(1.0), VecN(scale));
    return (translationMatrix * rotationMatrix * scaleMatrix);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
MatXxX<N + 1, RealType> Animation<N, RealType>::getInvTransformation(UInt frame, RealType frameFraction /*= RealType(0)*/)
{
    __NT_REQUIRE(m_bReady);
    return glm::inverse(getTransformation(frame, frameFraction));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
RealType Animation<N, RealType>::getUniformScale(UInt frame, RealType frameFraction /*= RealType(0)*/)
{
    if(m_KeyFrames.size() == 1 || (m_bPeriodic && frame < m_StartFrame)) {
        return m_KeyFrames[0].uniformScale;
    }

    if(m_bPeriodic && frame >= m_EndFrame) {
        frame = ((frame - m_StartFrame) % (m_EndFrame - m_StartFrame)) + m_StartFrame;
    }
    RealType x = static_cast<RealType>(frame) + frameFraction;
    return m_ScaleInterpolator(x);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(Animation)
