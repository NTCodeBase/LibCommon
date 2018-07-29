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

#include <Utils/NumberHelpers.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NumberHelpers
{
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define __BNN_INSTANTIATE_IRAND(type)                \
    template class MT_iRandom<type>;                 \
    template class iRand<type>;                      \
    template auto iRand<type>::vrnd<Vec2<type>>();   \
    template auto iRand<type>::vrnd<Vec3<type>>();   \
    template auto iRand<type>::mrnd<Mat2x2<type>>(); \
    template auto iRand<type>::mrnd<Mat3x3<type>>();

__BNN_INSTANTIATE_IRAND(Int)
__BNN_INSTANTIATE_IRAND(UInt)

template class MT_fRandom<float>;
template class MT_fRandom<double>;

#define __BNN_INSTANTIATE_FRAND(ClassName, type)         \
    template class ClassName<type>;                      \
    template auto ClassName<type>::vrnd<Vec2<type>>();   \
    template auto ClassName<type>::vrnd<Vec3<type>>();   \
    template auto ClassName<type>::mrnd<Mat2x2<type>>(); \
    template auto ClassName<type>::mrnd<Mat3x3<type>>();

__BNN_INSTANTIATE_FRAND(fRand,   float)
__BNN_INSTANTIATE_FRAND(fRand01, float)
__BNN_INSTANTIATE_FRAND(fRand11, float)

__BNN_INSTANTIATE_FRAND(fRand,   double)
__BNN_INSTANTIATE_FRAND(fRand01, double)
__BNN_INSTANTIATE_FRAND(fRand11, double)
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
}   // end namespace NumberHelpers
