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
#include <LibCommon/Data/DataIO.h>
#include <LibCommon/Utils/FileHelpers.h>
#include <LibCommon/Utils/NumberHelpers.h>

#include <algorithm>
#include <cassert>
#include <vector>
#include <istream>
#include <sstream>
#include <cstdlib>

#include <morton.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NTCodeBase {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T, bool USE_Z_ORDER = false>
class Array final {
public:
    using iterator               = typename StdVT<T>::iterator;
    using const_iterator         = typename StdVT<T>::const_iterator;
    using reverse_iterator       = typename StdVT<T>::reverse_iterator;
    using const_reverse_iterator = typename StdVT<T>::const_reverse_iterator;
    ////////////////////////////////////////////////////////////////////////////////
    // constructors & destructor
    Array() = default;
    Array(const Array<N, T, USE_Z_ORDER>& other) : m_Resolution(other.m_Resolution), m_Data(other.m_Data) {}

    template<class IndexType>
    Array(const VecX<N, IndexType>& size) {
        if constexpr (USE_Z_ORDER) {
            const auto maxSize = glm::compMax(size);
            m_Resolution = VecX<N, size_t>(NumberHelpers::nextPowerOfTwo(static_cast<size_t>(maxSize)));
        } else {
            for(Int d = 0; d < N; ++d) {
                m_Resolution[d] = static_cast<size_t>(size[d]);
            }
        }
        m_Data.resize(glm::compMul(m_Resolution));
    }

    ////////////////////////////////////////////////////////////////////////////////
    // assignment operator
    Array<N, T, USE_Z_ORDER>& operator=(const Array<N, T, USE_Z_ORDER>& other) {
        // check for self-assignment
        if(&other == this) {
            return *this;
        }
        m_Resolution = other.m_Resolution;
        m_Data       = other.m_Data;
        return *this;
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Array2D constructor =>
    template<class IndexType>
    Array(IndexType sizeX, IndexType sizeY) : m_Resolution(sizeX, sizeY), m_Data(sizeX * sizeY) {
        static_assert(N == 2, "Array dimension != 2");
    }

    template<class IndexType>
    Array(IndexType sizeX, IndexType sizeY, const StdVT<T>& data) : m_Resolution(sizeX, sizeY), m_Data(data) {
        static_assert(N == 2, "Array dimension != 2");
    }

    template<class IndexType>
    Array(IndexType sizeX, IndexType sizeY, const T& value) : m_Resolution(sizeX, sizeY), m_Data(sizeX * sizeY, value) {
        static_assert(N == 2, "Array dimension != 2");
    }

    template<class IndexType>
    Array(IndexType sizeX, IndexType sizeY, T* copyData) : m_Resolution(sizeX, sizeY), m_Data(sizeX * sizeY, copyData) {
        static_assert(N == 2, "Array dimension != 2");
    }

    // => end Array2D constructor
    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////
    // Array3D constructor =>
    template<class IndexType>
    Array(IndexType sizeX, IndexType sizeY, IndexType sizeZ) : m_Resolution(sizeX, sizeY, sizeZ), m_Data(sizeX * sizeY * sizeZ) {
        static_assert(N == 3, "Array dimension != 3");
    }

    template<class IndexType>
    Array(IndexType sizeX, IndexType sizeY, IndexType sizeZ, StdVT<T>& data) : m_Resolution(sizeX, sizeY, sizeZ), m_Data(data) {
        static_assert(N == 3, "Array dimension != 3");
    }

    template<class IndexType>
    Array(IndexType sizeX, IndexType sizeY, IndexType sizeZ, const T& value) : m_Resolution(sizeX, sizeY, sizeZ), m_Data(sizeX * sizeY * sizeZ, value) {
        static_assert(N == 3, "Array dimension != 3");
    }

    template<class IndexType>
    Array(IndexType sizeX, IndexType sizeY, IndexType sizeZ, T* copyData) : m_Resolution(sizeX, sizeY, sizeZ), m_Data(sizeX * sizeY * sizeZ, copyData) {
        static_assert(N == 3, "Array dimension != 3");
    }

    // => end Array3D constructor
    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////
    // size access
    size_t        empty() const { return m_Data.empty(); }
    size_t        capacity(void) const { return m_Data.capacity(); }
    size_t        dataSize(void) const { return m_Data.size(); }
    const VecX<N, size_t>& resolution() const { return m_Resolution; }

    ////////////////////////////////////////////////////////////////////////////////
    // index processing
    template<class IndexType>
    bool isValidIndex(const VecX<N, IndexType>& index) const {
        if(index[0] < 0 || static_cast<size_t>(index[0]) >= m_Resolution[0] ||
           index[1] < 0 || static_cast<size_t>(index[1]) >= m_Resolution[1]) {
            return false;
        }
        if constexpr (N == 3) {
            return index[2] >= 0 && static_cast<size_t>(index[2]) < m_Resolution[2];
        } else {
            return true;
        }
    }

    template<class IndexType>
    void checkIndex(const VecX<N, IndexType>& index) const {
        bool bIndexValid = isValidIndex<IndexType>(index);
        if(!bIndexValid) {
            std::stringstream ss;
            ss << "Invalid array index: ";

            for(Int d = 0; d < N - 1; ++d) {
                ss << index[d] << "/" << m_Resolution[d] << ", ";
            }
            ss << index[N - 1] << "/" << m_Resolution[N - 1];

            printf("%s\n", ss.str().c_str());
        }
        NT_REQUIRE(bIndexValid)
    }

    template<class IndexType>
    void checkZOrderIndex(size_t z_index, IndexType i, IndexType j) const {
        static_assert(N == 2, "Array dimension != 2");
        bool bIndexValid = z_index < m_Data.size();
        if(!bIndexValid) {
            printf("Invalid z_index = %zu, x = %u, y = %u, total size = %zu * %zu (%zu)\n", z_index, static_cast<UInt>(i), static_cast<UInt>(j),
                   m_Resolution[0], m_Resolution[1], m_Data.size());
        }
        NT_REQUIRE(bIndexValid)
    }

    template<class IndexType>
    void checkZOrderIndex(size_t z_index, IndexType i, IndexType j, IndexType k) const {
        static_assert(N == 3, "Array dimension != 3");
        bool bIndexValid = z_index < m_Data.size();
        if(!bIndexValid) {
            printf("Invalid z_index = %zu, x = %u, y = %u, z = %u, total size = %zu * %zu * %zu (%zu)\n", z_index,
                   static_cast<UInt>(i), static_cast<UInt>(j), static_cast<UInt>(k), m_Resolution[0], m_Resolution[1], m_Resolution[2], m_Data.size());
        }
        NT_REQUIRE(bIndexValid)
    }

    bool equalSize(const Array<N, T, USE_Z_ORDER>& other) const {
        for(Int d = 0; d < N; ++d) {
            if(m_Resolution[d] != other.m_Resolution[d]) {
                return false;
            }
        }
        return true;
    }

    template<class IndexType>
    bool equalSize(const VecX<N, IndexType>& otherSize) const {
        for(Int d = 0; d < N; ++d) {
            if(m_Resolution[d] != otherSize[d]) {
                return false;
            }
        }
        return true;
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Array2D =>
    template<class IndexType>
    bool isValidIndex(IndexType i, IndexType j) const {
        static_assert(N == 2, "Array dimension != 2");
        return (i >= 0 &&
                j >= 0 &&
                static_cast<size_t>(i) < m_Resolution[0] &&
                static_cast<size_t>(j) < m_Resolution[1]);
    }

    template<class IndexType>
    size_t getFlatIndex(IndexType i, IndexType j) const {
        static_assert(N == 2, "Array dimension != 2");
#ifndef NDEBUG
        checkIndex<IndexType>(Vec2<IndexType>(i, j));
#endif
        if constexpr (USE_Z_ORDER) {
            const auto z_index = libmorton::morton2D_64_encode(static_cast<uint_fast32_t>(i), static_cast<uint_fast32_t>(j));
#ifndef NDEBUG
            checkZOrderIndex(z_index, i, j);
#endif
            return z_index;
        } else {
            return static_cast<size_t>(j) * m_Resolution[0] + static_cast<size_t>(i);
        }
    }

    template<class IndexType>
    size_t getFlatIndex(const Vec2<IndexType>& index) const {
        static_assert(N == 2, "Array dimension != 2");
        return getFlatIndex<IndexType>(index[0], index[1]);
    }

    template<class IndexType>
    const T& operator()(IndexType i, IndexType j) const {
        static_assert(N == 2, "Array dimension != 2");
        return m_Data[getFlatIndex<IndexType>(i, j)];
    }

    template<class IndexType>
    T& operator()(IndexType i, IndexType j) {
        static_assert(N == 2, "Array dimension != 2");
        return m_Data[getFlatIndex<IndexType>(i, j)];
    }

    template<class IndexType>
    const T& operator()(const Vec2<IndexType>& index) const {
        static_assert(N == 2, "Array dimension != 2");
        return m_Data[getFlatIndex<IndexType>(index[0], index[1])];
    }

    template<class IndexType>
    T& operator()(const Vec2<IndexType>& index) {
        static_assert(N == 2, "Array dimension != 2");
        return m_Data[getFlatIndex<IndexType>(index[0], index[1])];
    }

    // => Array2D
    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////
    // Array3D =>
    template<class IndexType>
    bool isValidIndex(IndexType i, IndexType j, IndexType k) const {
        static_assert(N == 3, "Array dimension != 3");
        return (i >= 0 &&
                j >= 0 &&
                k >= 0 &&
                static_cast<size_t>(i) < m_Resolution[0] &&
                static_cast<size_t>(j) < m_Resolution[1] &&
                static_cast<size_t>(k) < m_Resolution[2]);
    }

    template<class IndexType>
    size_t getFlatIndex(IndexType i, IndexType j, IndexType k) const {
        static_assert(N == 3, "Array dimension != 3");
#ifndef NDEBUG
        checkIndex<IndexType>(Vec3<IndexType>(i, j, k));
#endif
        if constexpr (USE_Z_ORDER) {
            const auto z_index = libmorton::morton3D_64_encode(static_cast<uint_fast32_t>(i), static_cast<uint_fast32_t>(j), static_cast<uint_fast32_t>(k));
#ifndef NDEBUG
            checkZOrderIndex(z_index, i, j, k);
#endif
            return z_index;
        } else {
            return (static_cast<size_t>(k) * m_Resolution[1] + static_cast<size_t>(j)) * m_Resolution[0] + static_cast<size_t>(i);
        }
    }

    template<class IndexType>
    size_t getFlatIndex(const Vec3<IndexType>& index) const {
        static_assert(N == 3, "Array dimension != 3");
        return getFlatIndex<IndexType>(index[0], index[1], index[2]);
    }

    template<class IndexType>
    const T& operator()(IndexType i, IndexType j, IndexType k) const {
        static_assert(N == 3, "Array dimension != 3");
        return m_Data[getFlatIndex<IndexType>(i, j, k)];
    }

    template<class IndexType>
    T& operator()(IndexType i, IndexType j, IndexType k) {
        static_assert(N == 3, "Array dimension != 3");
        return m_Data[getFlatIndex<IndexType>(i, j, k)];
    }

    template<class IndexType>
    const T& operator()(const Vec3<IndexType>& index) const {
        static_assert(N == 3, "Array dimension != 3");
        return m_Data[getFlatIndex<IndexType>(index[0], index[1], index[2])];
    }

    template<class IndexType>
    T& operator()(const Vec3<IndexType>& index) {
        static_assert(N == 3, "Array dimension != 3");
        return m_Data[getFlatIndex<IndexType>(index[0], index[1], index[2])];
    }

    // => Array3D
    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////
    // flat data
    StdVT<T>& flatData() { return m_Data; }
    const StdVT<T>& flatData() const { return m_Data; }

    template<class IndexType> T& flatData(IndexType flatIdx) { return m_Data[flatIdx]; }
    template<class IndexType> const T& flatData(IndexType flatIdx) const { return m_Data[flatIdx]; }
    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////
    // iterators
    const auto& front(void) const { assert(m_Data.size() > 0); return m_Data.front(); }
    auto& front(void) { assert(m_Data.size() > 0); return m_Data.front(); }

    const auto& back(void) const { assert(m_Data.size() > 0); return m_Data.back(); }
    auto& back(void) { assert(m_Data.size() > 0); return m_Data.back(); }

    auto begin(void) { return m_Data.begin(); }
    auto end(void) { return m_Data.end(); }

    auto cbegin(void) const { return m_Data.cbegin(); }
    auto cend(void) const { return m_Data.cend(); }

    auto rbegin(void) { return reverse_iterator(end()); }
    auto rend(void) { return reverse_iterator(begin()); }

    auto crbegin(void) const { return const_reverse_iterator(cend()); }
    auto crend(void) const { return const_reverse_iterator(cbegin()); }

    ////////////////////////////////////////////////////////////////////////////////
    // data manipulation
    template<class IndexType>
    void assign(const VecX<N, IndexType>& size, const T& value) {
        m_Resolution = size;
        m_Data.assign(glm::compMul(m_Resolution), value);
    }

    template<class IndexType>
    void assign(const VecX<N, IndexType>& size, const T* copydata) {
        m_Resolution = size;
        m_Data.assign(glm::compMul(m_Resolution), copydata);
    }

    void assign(const T& value) { m_Data.assign(m_Data.size(), value); }
    void copyDataFrom(const Array<N, T, USE_Z_ORDER>& other) { NT_REQUIRE(equalSize(other)) m_Data = other.m_Data; }
    void setZero() { m_Data.assign(m_Data.size(), 0); }
    void clear() { m_Data.resize(0); m_Resolution = VecX<N, size_t>(0); }
    void swap(Array<N, T, USE_Z_ORDER>& other) { std::swap(m_Resolution, other.m_Resolution); m_Data.swap(other.m_Data); }

    template<class IndexType>
    void reserve(IndexType size) { m_Data.reserve(size); }

    template<class IndexType>
    void reserve(const VecX<N, IndexType>& size) {
        if constexpr (USE_Z_ORDER) {
            auto      maxSize = NumberHelpers::nextPowerOfTwo(glm::compMax(size));
            IndexType size1D  = IndexType(1);
            for(int i = 0; i < N; ++i) {
                size1D *= maxSize;
            }
            m_Data.reserve(size1D);
        } else {
            m_Data.reserve(glm::compMul(size));
        }
    }

    template<class IndexType>
    void resize(const VecX<N, IndexType>& newSize) {
        if constexpr (USE_Z_ORDER) {
            auto maxSize = glm::compMax(newSize);
            m_Resolution = VecX<N, size_t>(NumberHelpers::nextPowerOfTwo(static_cast<size_t>(maxSize)));
        } else {
            m_Resolution = newSize;
        }
        m_Data.resize(glm::compMul(m_Resolution));
    }

    template<class IndexType>
    void resize(const VecX<N, IndexType>& newSize, const T& value) {
        if constexpr (USE_Z_ORDER) {
            auto maxSize = glm::compMax(newSize);
            m_Resolution = VecX<N, size_t>(NumberHelpers::nextPowerOfTwo(static_cast<size_t>(maxSize)));
        } else {
            m_Resolution = newSize;
        }
        m_Data.resize(glm::compMul(m_Resolution), value);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Array2D =>
    template<class IndexType>
    void assign(IndexType sizeX, IndexType sizeY, const T& value) {
        static_assert(N == 2, "Array dimension != 2");
        assign(Vec2<IndexType>(sizeX, sizeY), value);
    }

    template<class IndexType>
    void assign(IndexType sizeX, IndexType sizeY, const T* copydata) {
        static_assert(N == 2, "Array dimension != 2");
        assign(Vec2<IndexType>(sizeX, sizeY), copydata);
    }

    template<class IndexType>
    void reserve(IndexType sizeX, IndexType sizeY) {
        static_assert(N == 2, "Array dimension != 2");
        reserve(Vec2<IndexType>(sizeX, sizeY));
    }

    template<class IndexType>
    void resize(IndexType sizeX, IndexType sizeY) {
        static_assert(N == 2, "Array dimension != 2");
        resize(Vec2<IndexType>(sizeX, sizeY));
    }

    template<class IndexType>
    void resize(IndexType sizeX, IndexType sizeY, const T& value) {
        static_assert(N == 2, "Array dimension != 2");
        resize(Vec2<IndexType>(sizeX, sizeY), value);
    }

    // => Array2D
    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////
    // Array3D =>
    template<class IndexType>
    void assign(IndexType sizeX, IndexType sizeY, IndexType sizeZ, const T& value) {
        static_assert(N == 3, "Array dimension != 3");
        assign(Vec3<IndexType>(sizeX, sizeY, sizeZ), value);
    }

    template<class IndexType>
    void assign(IndexType sizeX, IndexType sizeY, IndexType sizeZ, const T* copydata) {
        static_assert(N == 3, "Array dimension != 3");
        assign(Vec3<IndexType>(sizeX, sizeY, sizeZ), copydata);
    }

    template<class IndexType>
    void reserve(IndexType sizeX, IndexType sizeY, IndexType sizeZ) {
        static_assert(N == 3, "Array dimension != 3");
        reserve(Vec2<IndexType>(sizeX, sizeY, sizeZ));
    }

    template<class IndexType>
    void resize(IndexType sizeX, IndexType sizeY, IndexType sizeZ) {
        static_assert(N == 3, "Array dimension != 3");
        resize(Vec3<IndexType>(sizeX, sizeY, sizeZ));
    }

    template<class IndexType>
    void resize(IndexType sizeX, IndexType sizeY, IndexType sizeZ, const T& value) {
        static_assert(N == 3, "Array dimension != 3");
        resize(Vec3<IndexType>(sizeX, sizeY, sizeZ), value);
    }

    // => Array2D
    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////
    // file IO
    bool saveToFile(const String& fileName) {
        DataBuffer buffer;
        for(Int d = 0; d < N; ++d) {
            buffer.append<UInt>(static_cast<UInt>(m_Resolution[d]));
        }
        buffer.append(m_Data, false);
        return FileHelpers::writeFile(buffer.data(), buffer.size(), fileName);
    }

    bool loadFromFile(const String& fileName) {
        DataBuffer buffer;
        if(!FileHelpers::readFile(buffer.buffer(), fileName)) {
            return false;
        }

        for(UInt d = 0; d < N; ++d) {
            UInt tmp;
            buffer.getData<UInt>(tmp, sizeof(UInt) * d);
            m_Resolution[d] = static_cast<size_t>(tmp);
        }
        NT_REQUIRE(buffer.buffer().size() == N * sizeof(UInt) + sizeof(T) * glm::compMul(m_Resolution))
        buffer.getData(m_Data, sizeof(UInt) * N, static_cast<UInt>(glm::compMul(m_Resolution)));
        return true;
    }

private:
    VecX<N, size_t> m_Resolution = VecX<N, size_t>(0);
    StdVT<T>        m_Data;
}; // end class Array

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T> using Array2 = Array<2, T>;
template<class T> using Array3 = Array<3, T>;
////////////////////////////////////////////////////////////////////////////////
using Array2c  = Array2<char>;
using Array2uc = Array2<unsigned char>;
using Array2s  = Array2<short>;
using Array2i  = Array2<int>;
using Array2ui = Array2<unsigned int>;
using Array2f  = Array2<float>;
using Array2d  = Array2<double>;
////////////////////////////////////////////////////////////////////////////////
using Array3c  = Array3<char>;
using Array3uc = Array3<unsigned char>;
using Array3s  = Array3<short>;
using Array3us = Array3<unsigned short>;
using Array3i  = Array3<int>;
using Array3ui = Array3<unsigned int>;
using Array3f  = Array3<float>;
using Array3d  = Array3<double>;
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace NTCodeBase
