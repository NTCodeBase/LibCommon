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

#include <LibCommon/Animation/Animation.h>
#include <algorithm>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void RigidBodyAnimation<N, Real_t>::makeReady(bool bCubicIntTranslation, bool bCubicIntRotation)
{
    if(m_KeyFrames.size() == 0) {
        return;
    }
    ////////////////////////////////////////////////////////////////////////////////
    // sort key frames by frame order
    std::sort(m_KeyFrames.begin(), m_KeyFrames.end(), [](const auto& k1, const auto& k2) { return k1.frame < k2.frame; });
    ////////////////////////////////////////////////////////////////////////////////
    // add more 1 point if there are only 2 key frames
    if(m_KeyFrames.size() == 2) {
        auto& keyFrame  = m_KeyFrames.emplace_back(KeyFrame<N, Real_t>());
        auto& keyFrame0 = m_KeyFrames[0];
        auto& keyFrame1 = m_KeyFrames[1];
        keyFrame.frame        = (keyFrame0.frame + keyFrame1.frame) / 2;
        keyFrame.translation  = (keyFrame0.translation + keyFrame1.translation) * Real_t(0.5);
        keyFrame.rotation     = (keyFrame0.rotation + keyFrame1.rotation) * Real_t(0.5);
        keyFrame.uniformScale = (keyFrame0.uniformScale + keyFrame1.uniformScale) * Real_t(0.5);
        keyFrame.invScale     = Real_t(1) / keyFrame.uniformScale;
        std::swap(keyFrame, keyFrame1);
    }
    ////////////////////////////////////////////////////////////////////////////////
    m_StartFrame = m_KeyFrames.front().frame;
    m_EndFrame   = m_KeyFrames.back().frame;
    __NT_REQUIRE(m_EndFrame > m_StartFrame);
    m_FrameRange = m_EndFrame - m_StartFrame;
    ////////////////////////////////////////////////////////////////////////////////
    StdVT<Real_t> frames;
    StdVT<Real_t> translations[N];
    StdVT<Real_t> rotations[N + 1];

    frames.reserve(nKeyFrames());
    for(Int d = 0; d < N; ++d) {
        translations[d].reserve(nKeyFrames());
        rotations[d].reserve(nKeyFrames());
    }
    rotations[N].reserve(nKeyFrames());
    ////////////////////////////////////////////////////////////////////////////////
    for(const auto& keyFrame : m_KeyFrames) {
        frames.push_back(static_cast<Real_t>(keyFrame.frame));
        for(Int d = 0; d < N; ++d) {
            translations[d].push_back(keyFrame.translation[d]);
            rotations[d].push_back(keyFrame.rotation[d]);
        }
        rotations[N].push_back(keyFrame.rotation[N]);
    }
    ////////////////////////////////////////////////////////////////////////////////
    for(Int d = 0; d < N; ++d) {
        m_TranslationInterpolator[d].setBoundary(CubicSpline<Real_t>::BDType::FirstOrder, 0, CubicSpline<Real_t>::BDType::FirstOrder, 0);
        m_TranslationInterpolator[d].setPoints(frames, translations[d], bCubicIntTranslation);

        m_RotationInterpolator[d].setBoundary(CubicSpline<Real_t>::BDType::FirstOrder, 0, CubicSpline<Real_t>::BDType::FirstOrder, 0);
        m_RotationInterpolator[d].setPoints(frames, rotations[d], bCubicIntRotation);
    }
    m_RotationInterpolator[N].setBoundary(CubicSpline<Real_t>::BDType::FirstOrder, 0, CubicSpline<Real_t>::BDType::FirstOrder, 0);
    m_RotationInterpolator[N].setPoints(frames, rotations[N], bCubicIntRotation);
    ////////////////////////////////////////////////////////////////////////////////
    m_bReady = true;
    m_StartFrameTransformationMatrix = getTransformationMatrix(m_StartFrame);
    m_EndFrameTransformationMatrix   = getTransformationMatrix(m_EndFrame);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
MatXxX<N + 1, Real_t> RigidBodyAnimation<N, Real_t>::getInactiveTransformationMatrix(UInt frame)
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
template<Int N, class Real_t>
MatXxX<N + 1, Real_t> RigidBodyAnimation<N, Real_t>::getTransformationMatrix(UInt frame, Real_t frameFraction /*= Real_t(0)*/)
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
            frameFraction = Real_t(0);
        }
    }

    Real_t x = static_cast<Real_t>(frame) + frameFraction;
    for(Int d = 0; d < N; ++d) {
        translation[d] = m_TranslationInterpolator[d](x);
        rotation[d]    = m_RotationInterpolator[d](x);
    }
    rotation[N] = m_RotationInterpolator[N](x);

    MatNp1xNp1 transMatrix(1);
    if(rotation[N] != 0) {
        transMatrix = glm::rotate(transMatrix, rotation[N], VecN(rotation));
    }
    transMatrix = glm::translate(transMatrix, translation);
    return transMatrix;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
MatXxX<N + 1, Real_t> RigidBodyAnimation<N, Real_t>::getInvTransformation(UInt frame, Real_t frameFraction /*= Real_t(0)*/)
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
template<Int N, class Real_t>
void Animation<N, Real_t>::makeReady(bool bCubicIntTranslation, bool bCubicIntRotation, bool bCubicIntScale)
{
    if(this->m_KeyFrames.size() == 0) {
        return;
    }
    ////////////////////////////////////////////////////////////////////////////////
    StdVT<Real_t> frames;
    StdVT<Real_t> translations[N];
    StdVT<Real_t> rotations[N + 1];
    StdVT<Real_t> scales;

    frames.reserve(this->nKeyFrames());
    for(Int d = 0; d < N; ++d) {
        translations[d].reserve(this->nKeyFrames());
        rotations[d].reserve(this->nKeyFrames());
    }
    rotations[N].reserve(this->nKeyFrames());
    scales.reserve(this->nKeyFrames());
    ////////////////////////////////////////////////////////////////////////////////
    UInt maxFrame = 0;
    for(const auto& keyFrame : this->m_KeyFrames) {
        frames.push_back(static_cast<Real_t>(keyFrame.frame));
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
    if(this->m_EndFrame == 0) {
        // if end frame has not been set, set it to the latest key frame
        this->m_EndFrame = maxFrame;
    }
    ////////////////////////////////////////////////////////////////////////////////
    for(Int d = 0; d < N; ++d) {
        this->m_TranslationInterpolator[d].setBoundary(CubicSpline<Real_t>::BDType::FirstOrder, 0, CubicSpline<Real_t>::BDType::FirstOrder, 0);
        this->m_TranslationInterpolator[d].setPoints(frames, translations[d], bCubicIntTranslation);

        this->m_RotationInterpolator[d].setBoundary(CubicSpline<Real_t>::BDType::FirstOrder, 0, CubicSpline<Real_t>::BDType::FirstOrder, 0);
        this->m_RotationInterpolator[d].setPoints(frames, rotations[d], bCubicIntRotation);
    }
    this->m_RotationInterpolator[N].setBoundary(CubicSpline<Real_t>::BDType::FirstOrder, 0, CubicSpline<Real_t>::BDType::FirstOrder, 0);
    this->m_RotationInterpolator[N].setPoints(frames, rotations[N], bCubicIntRotation);
    this->m_ScaleInterpolator.setBoundary(CubicSpline<Real_t>::BDType::FirstOrder, 0, CubicSpline<Real_t>::BDType::FirstOrder, 0);
    this->m_ScaleInterpolator.setPoints(frames, scales, bCubicIntScale);
    ////////////////////////////////////////////////////////////////////////////////
    __NT_REQUIRE(this->m_EndFrame > this->m_StartFrame);
    this->m_FrameRange = this->m_EndFrame - this->m_StartFrame;
    this->m_bReady     = true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
Real_t Animation<N, Real_t>::getUniformScale(UInt frame, Real_t frameFraction /*= Real_t(0)*/)
{
    if(this->nKeyFrames() == 0) {
        return Real_t(1);
    }
    ////////////////////////////////////////////////////////////////////////////////
    if(frame < this->m_StartFrame) {
        frame = this->m_StartFrame;
    } else if(frame >= this->m_EndFrame) {
        if(this->m_bPeriodic) {
            frame = ((frame - this->m_StartFrame) % this->m_FrameRange) + this->m_StartFrame;
        } else {
            frame         = this->m_EndFrame;
            frameFraction = Real_t(0);
        }
    }
    Real_t x = static_cast<Real_t>(frame) + frameFraction;
    return m_ScaleInterpolator(x);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
MatXxX<N + 1, Real_t> Animation<N, Real_t>::getTransformationMatrix(UInt frame, Real_t frameFraction /*= Real_t(0)*/)
{
    if(this->nKeyFrames() == 0) {
        return MatNp1xNp1(1);
    }
    ////////////////////////////////////////////////////////////////////////////////
    __NT_REQUIRE(this->m_bReady)
    VecN translation;
    VecNp1 rotation;
    Real_t scale;

    if(frame < this->m_StartFrame) {
        frame = this->m_StartFrame;
    } else if(frame >= this->m_EndFrame) {
        if(this->m_bPeriodic) {
            frame = ((frame - this->m_StartFrame) % this->m_FrameRange) + this->m_StartFrame;
        } else {
            frame         = this->m_EndFrame;
            frameFraction = Real_t(0);
        }
    }

    Real_t x = static_cast<Real_t>(frame) + frameFraction;
    for(Int d = 0; d < N; ++d) {
        translation[d] = this->m_TranslationInterpolator[d](x);
        rotation[d]    = this->m_RotationInterpolator[d](x);
    }
    rotation[N] = this->m_RotationInterpolator[N](x);
    scale       = this->m_ScaleInterpolator(x);

    MatNp1xNp1 transMatrix(1);
    transMatrix = glm::scale(transMatrix, VecN(scale));
    if(rotation[N] != 0) {
        transMatrix = glm::rotate(transMatrix, rotation[N], VecN(rotation));
    }
    transMatrix = glm::translate(transMatrix, translation);
    return transMatrix;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(RigidBodyAnimation)
