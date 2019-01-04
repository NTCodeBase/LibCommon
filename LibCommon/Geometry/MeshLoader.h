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

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NTCodeBase {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// This class need to be a template, as it can be used in simulation(double) and rendering(float)
class MeshLoader {
public:
    MeshLoader() { clearData(); }
    MeshLoader(String meshFile) { loadMesh(meshFile); }

    bool loadMesh(const String& meshFile);
    bool isMeshReady() const { return m_isMeshReady; }
    void scaleToBox();

    auto  getMeshCenter() const { assert(m_isMeshReady); return float(0.5) * (m_AABBMin + m_AABBMax); }
    auto  getNTriangles() const { return m_NumTriangles; }
    const auto& getAABBMin() const { assert(m_isMeshReady); return m_AABBMin; }
    const auto& getAABBMax() const { assert(m_isMeshReady); return m_AABBMax; }

    Vec3f getCameraPosition(Vec3f camDirection, float fov = 45);
    float getCameraDistance(float fov);

    const auto& getVertices() const { assert(m_isMeshReady); return m_Vertices; }
    const auto& getNormals() const { assert(m_isMeshReady); return m_Normals; }
    const auto& getTexCoord2D() const { assert(m_isMeshReady); return m_TexCoord2D; }

    const auto& getFaces() const { assert(m_isMeshReady); return m_Faces; }
    const auto& getFaceVertices() const { assert(m_isMeshReady); return m_FaceVertices; }
    const auto& getFaceVertexNormals() const { assert(m_isMeshReady); return m_FaceVertexNormals; }
    const auto& getFaceVertexColors() const { assert(m_isMeshReady); return m_FaceVertexColors; }
    const auto& getFaceVTexCoord2D() const { assert(m_isMeshReady); return m_FaceVertexTexCoord2D; }

    auto getNFaces() const noexcept { assert(m_isMeshReady); return (m_Faces.size() / 3); }
    auto getNVertices() const noexcept { assert(m_isMeshReady); return (m_Vertices.size() / 3); }
    auto getNFaceVertices() const noexcept { assert(m_isMeshReady); return (m_FaceVertices.size() / 3); }

    void swapXY() { swapCoordinates(0, 1); }
    void swapYZ() { swapCoordinates(1, 2); }
    void swapXZ() { swapCoordinates(0, 2); }
private:
    enum class MeshFileType {
        OBJFile,
        UnsupportedFileType
    };
    MeshFileType getMeshType(const String& meshFile);
    bool         loadObj(const String& meshFile);

    void clearData();
    void swapCoordinates(int k1, int k2);
    void computeFaceVertexData();
    ////////////////////////////////////////////////////////////////////////////////

    UInt m_NumTriangles = 0;
    bool m_isMeshReady  = false;

    Vec3f m_AABBMin;
    Vec3f m_AABBMax;

    StdVT_UInt  m_Faces;
    StdVT_Float m_Vertices;
    StdVT_Float m_Normals;
    StdVT_Float m_TexCoord2D;
    StdVT_Float m_FaceVertices;
    StdVT_Float m_FaceVertexNormals;
    StdVT_Float m_FaceVertexColors;
    StdVT_Float m_FaceVertexTexCoord2D;
};
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace NTCodeBase
