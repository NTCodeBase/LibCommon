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
void Animation<N, RealType>::addKeyFrame(UInt frame, const VecX<N + 1, RealType>& rotation)
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
void Animation<N, RealType>::addKeyFrame(UInt frame, const VecN& translation, const VecX<N + 1, RealType>& rotation, RealType scale)
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
    size_t nKeyFrames = m_KeyFrames.size();

    StdVT<RealType> frames;
    StdVT<RealType> translations[N];
    StdVT<RealType> rotations[N + 1];
    StdVT<RealType> scales;

    for(Int d = 0; d < N; ++d) {
        translations[d].reserve(nKeyFrames);
        rotations[d].reserve(nKeyFrames);
    }
    rotations[N].reserve(nKeyFrames);
    scales.reserve(nKeyFrames);
    frames.reserve(nKeyFrames);

    ////////////////////////////////////////////////////////////////////////////////
    for(const auto& keyFrame : m_KeyFrames) {
        if(m_bPeriodic && keyFrame.frame < m_StartFrame) {
            continue;
        }
        m_MaxFrame = (m_MaxFrame < keyFrame.frame) ? keyFrame.frame : m_MaxFrame;

        for(Int d = 0; d < N; ++d) {
            translations[d].push_back(keyFrame.translation[d]);
            rotations[d].push_back(keyFrame.rotation[d]);
        }
        rotations[N].push_back(keyFrame.rotation[N]);
        scales.push_back(keyFrame.uniformScale);
        frames.push_back(static_cast<RealType>(keyFrame.frame));
    }

    for(Int d = 0; d < N; ++d) {
        m_TranslationSpline[d].setBoundary(CubicSpline<RealType>::BDType::FirstOrder, 0, CubicSpline<RealType>::BDType::FirstOrder, 0);
        m_TranslationSpline[d].setPoints(frames, translations[d], bCubicIntTranslation);

        m_RotationSpline[d].setBoundary(CubicSpline<RealType>::BDType::FirstOrder, 0, CubicSpline<RealType>::BDType::FirstOrder, 0);
        m_RotationSpline[d].setPoints(frames, rotations[d], bCubicIntRotation);
    }
    m_RotationSpline[N].setBoundary(CubicSpline<RealType>::BDType::FirstOrder, 0, CubicSpline<RealType>::BDType::FirstOrder, 0);
    m_RotationSpline[N].setPoints(frames, rotations[N], bCubicIntRotation);
    m_ScaleSpline.setBoundary(CubicSpline<RealType>::BDType::FirstOrder, 0, CubicSpline<RealType>::BDType::FirstOrder, 0);
    m_ScaleSpline.setPoints(frames, scales, bCubicIntScale);

    ////////////////////////////////////////////////////////////////////////////////
    m_bReady = true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
void Animation<N, RealType>::getTransformation(VecN& translation, VecX<N + 1, RealType>& rotation, RealType& scale, UInt frame, RealType fraction)
{
    if(m_KeyFrames.size() == 1 || (m_bPeriodic && frame < m_StartFrame)) {
        translation = m_KeyFrames[0].translation;
        rotation    = m_KeyFrames[0].rotation;
        scale       = m_KeyFrames[0].uniformScale;

        return;
    }

    ////////////////////////////////////////////////////////////////////////////////
    __NT_REQUIRE(m_bReady)

    if(m_bPeriodic && frame > m_MaxFrame) {
        frame = ((frame - m_StartFrame) % (m_MaxFrame - m_StartFrame)) + m_StartFrame;
    }
    RealType x = static_cast<RealType>(frame) + fraction;

    for(Int d = 0; d < N; ++d) {
        translation[d] = m_TranslationSpline[d](x);
        rotation[d]    = m_RotationSpline[d](x);
    }
    rotation[N] = m_RotationSpline[N](x);
    scale       = m_ScaleSpline(x);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
MatXxX<N + 1, RealType> Animation<N, RealType>::getTransformation(UInt frame, RealType fraction)
{
    if(m_KeyFrames.size() == 1 || (m_bPeriodic && frame < m_StartFrame)) {
        MatXxX<N + 1, RealType> translationMatrix = glm::translate(MatXxX<N + 1, RealType>(1.0), m_KeyFrames[0].translation);
        MatXxX<N + 1, RealType> rotationMatrix    = glm::rotate(MatXxX<N + 1, RealType>(1.0), m_KeyFrames[0].rotation[N], VecN(m_KeyFrames[0].rotation));
        MatXxX<N + 1, RealType> scaleMatrix       = glm::scale(MatXxX<N + 1, RealType>(1.0), VecN(m_KeyFrames[0].uniformScale));

        return (translationMatrix * rotationMatrix * scaleMatrix);
    }

    ////////////////////////////////////////////////////////////////////////////////
    __NT_REQUIRE(m_bReady)
    VecN translation;
    VecX<N + 1, RealType> rotation;
    RealType              scale;

    if(m_bPeriodic && frame > m_MaxFrame) {
        frame = ((frame - m_StartFrame) % (m_MaxFrame - m_StartFrame)) + m_StartFrame;
    }
    RealType x = static_cast<RealType>(frame) + fraction;

    for(Int d = 0; d < N; ++d) {
        translation[d] = m_TranslationSpline[d](x);
        rotation[d]    = m_RotationSpline[d](x);
    }
    rotation[N] = m_RotationSpline[N](x);
    scale       = m_ScaleSpline(x);

    MatXxX<N + 1, RealType> translationMatrix = glm::translate(MatXxX<N + 1, RealType>(1.0), translation);
    MatXxX<N + 1, RealType> rotationMatrix    = glm::rotate(MatXxX<N + 1, RealType>(1.0), rotation[N], VecN(rotation));
    MatXxX<N + 1, RealType> scaleMatrix       = glm::scale(MatXxX<N + 1, RealType>(1.0), VecN(scale));

    return (translationMatrix * rotationMatrix * scaleMatrix);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
MatXxX<N + 1, RealType> Animation<N, RealType>::getInvTransformation(UInt frame, RealType fraction)
{
    __NT_REQUIRE(m_bReady);
    return glm::inverse(getTransformation(frame, fraction));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType>
RealType Animation<N, RealType>::getUniformScale(UInt frame, RealType fraction)
{
    if(m_KeyFrames.size() == 1) {
        return m_KeyFrames[0].uniformScale;
    }

    if(m_bPeriodic && frame > m_MaxFrame) {
        frame = frame % m_MaxFrame;
    }
    RealType x = static_cast<RealType>(frame) + fraction;
    return m_ScaleSpline(x);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(Animation)
