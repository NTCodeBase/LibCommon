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

#include <Geometry/MeshLoader.h>
#include <Utils/MathHelpers.h>

#include <tiny_obj_loader.h>

#include <fstream>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
bool MeshLoader::loadMesh(const String& meshFile)
{
    auto meshType = getMeshType(meshFile);
    if(meshType == MeshFileType::UnsupportedFileType) {
        return false;
    }

    ////////////////////////////////////////////////////////////////////////////////
    // => clear data
    clearData();

    ////////////////////////////////////////////////////////////////////////////////
    // => load mesh
    bool result = loadObj(meshFile);
    m_isMeshReady = result;
    if(result) {
        computeFaceVertexData();
    }
    return result;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void MeshLoader::scaleToBox()
{
    Vec3f diff    = m_AABBMax - m_AABBMin;
    float maxSize = fmaxf(fmaxf(fabsf(diff[0]), fabsf(diff[1])), fabsf(diff[2]));
    float scale   = 2.0f / maxSize;

    // multiply all vertices by scale to make the mesh having max(w, h, d) = 1
    m_AABBMin = m_AABBMin * scale;
    m_AABBMax = m_AABBMax * scale;

    // expand the bounding box
    Vec3f meshCenter = (m_AABBMax + m_AABBMin) * 0.5f;
    auto  cmin       = m_AABBMin - meshCenter;
    auto  cmax       = m_AABBMax - meshCenter;

    m_AABBMin = meshCenter + glm::normalize(cmin) * glm::length(cmin);
    m_AABBMax = meshCenter + glm::normalize(cmax) * glm::length(cmax);

    // to move the mesh center to origin
    m_AABBMin = m_AABBMin - meshCenter;
    m_AABBMax = m_AABBMax - meshCenter;

    Vec3f* vertexPtr = reinterpret_cast<Vec3f*>(m_Vertices.data());
    for(size_t i = 0, iend = m_Vertices.size() / 3; i < iend; ++i) {
        vertexPtr[i] = vertexPtr[i] * scale;
        vertexPtr[i] = vertexPtr[i] - meshCenter;
    }

    Vec3f* faceVertexPtr = reinterpret_cast<Vec3f*>(m_FaceVertices.data());
    for(size_t i = 0, iend = m_FaceVertices.size() / 3; i < iend; ++i) {
        faceVertexPtr[i] = faceVertexPtr[i] * scale;
        faceVertexPtr[i] = faceVertexPtr[i] - meshCenter;
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
Vec3f MeshLoader::getCameraPosition(Vec3f camDirection, float fov /*= 45*/)
{
    return camDirection * getCameraDistance(fov * float(0.75)) + getMeshCenter();
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
float MeshLoader::getCameraDistance(float fov)
{
    float halfLength = (m_AABBMax.y - m_AABBMin.y) * float(0.5);

    return (halfLength / std::tan(fov * float(0.5 * M_PI / 180.0)));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void MeshLoader::swapCoordinates(int k1, int k2)
{
    __NT_REQUIRE(k1 >= 0 && k1 <= 2 && k2 >= 0 && k2 <= 2 && k1 != k2);
    for(size_t i = 0, iend = m_Vertices.size() / 3; i < iend; ++i) {
        std::swap(m_Vertices[i * 3 + k1], m_Vertices[i * 3 + k2]);
    }

    for(size_t i = 0, iend = m_FaceVertices.size() / 3; i < iend; ++i) {
        std::swap(m_FaceVertices[i * 3 + k1], m_FaceVertices[i * 3 + k2]);
    }

    std::swap(m_AABBMin[k1], m_AABBMin[k2]);
    std::swap(m_AABBMax[k1], m_AABBMax[k2]);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
MeshLoader::MeshFileType MeshLoader::getMeshType(const String& meshFile)
{
    const std::string extension = meshFile.substr(meshFile.rfind('.') + 1);
    if(extension == "obj" || extension == "OBJ" || extension == "Obj") {
        return MeshFileType::OBJFile;
    }
    return MeshFileType::UnsupportedFileType;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void MeshLoader::clearData()
{
    m_isMeshReady  = false;
    m_NumTriangles = 0;

    m_AABBMin = Vec3f(1e10);
    m_AABBMax = Vec3f(-1e10);

    m_Vertices.resize(0);
    m_Normals.resize(0);

    m_Faces.resize(0);
    m_FaceVertices.resize(0);
    m_FaceVertexNormals.resize(0);
    m_FaceVertexColors.resize(0);
    m_FaceVertexTexCoord2D.resize(0);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
bool MeshLoader::loadObj(const String& meshFile)
{
    std::vector<tinyobj::shape_t>    obj_shapes;
    std::vector<tinyobj::material_t> obj_materials;
    tinyobj::attrib_t                attrib;

    std::string errStr;
    bool        result = tinyobj::LoadObj(&attrib, &obj_shapes, &obj_materials, &errStr, meshFile.c_str(), NULL, true);
    if(!errStr.empty()) {
        std::cerr << "tinyobj: " << errStr << std::endl;
    }

    if(!result) {
        std::cerr << "Failed to load " << meshFile << std::endl;
        return false;
    }

    m_Vertices   = attrib.vertices;
    m_Normals    = attrib.normals;
    m_TexCoord2D = attrib.texcoords;
    ////////////////////////////////////////////////////////////////////////////////
    // => convert data
    for(size_t s = 0; s < obj_shapes.size(); s++) {
        for(size_t f = 0; f < obj_shapes[s].mesh.indices.size() / 3; ++f) {
            ++m_NumTriangles;

            tinyobj::index_t idx0 = obj_shapes[s].mesh.indices[3 * f + 0];
            tinyobj::index_t idx1 = obj_shapes[s].mesh.indices[3 * f + 1];
            tinyobj::index_t idx2 = obj_shapes[s].mesh.indices[3 * f + 2];

            int v0 = idx0.vertex_index;
            int v1 = idx1.vertex_index;
            int v2 = idx2.vertex_index;
            assert(v0 >= 0);
            assert(v1 >= 0);
            assert(v2 >= 0);

            m_Faces.push_back(static_cast<UInt>(v0));
            m_Faces.push_back(static_cast<UInt>(v1));
            m_Faces.push_back(static_cast<UInt>(v2));

            Vec3f v[3];
            for(int k = 0; k < 3; ++k) {
                v[0][k] = m_Vertices[3 * v0 + k];
                v[1][k] = m_Vertices[3 * v1 + k];
                v[2][k] = m_Vertices[3 * v2 + k];

                m_AABBMin[k] = MathHelpers::min(v[0][k], m_AABBMin[k]);
                m_AABBMin[k] = MathHelpers::min(v[1][k], m_AABBMin[k]);
                m_AABBMin[k] = MathHelpers::min(v[2][k], m_AABBMin[k]);

                m_AABBMax[k] = MathHelpers::max(v[0][k], m_AABBMax[k]);
                m_AABBMax[k] = MathHelpers::max(v[1][k], m_AABBMax[k]);
                m_AABBMax[k] = MathHelpers::max(v[2][k], m_AABBMax[k]);
            }

            for(int k = 0; k < 3; ++k) {
                m_FaceVertices.push_back(v[k][0]);
                m_FaceVertices.push_back(v[k][1]);
                m_FaceVertices.push_back(v[k][2]);
            }

            if(attrib.normals.size() > 0) {
                Vec3f n[3];
                int   n0 = idx0.normal_index;
                int   n1 = idx1.normal_index;
                int   n2 = idx2.normal_index;
                assert(n0 >= 0);
                assert(n1 >= 0);
                assert(n2 >= 0);

                for(int k = 0; k < 3; ++k) {
                    n[0][k] = attrib.normals[3 * n0 + k];
                    n[1][k] = attrib.normals[3 * n1 + k];
                    n[2][k] = attrib.normals[3 * n2 + k];
                }

                for(int k = 0; k < 3; ++k) {
                    m_FaceVertexNormals.push_back(n[k][0]);
                    m_FaceVertexNormals.push_back(n[k][1]);
                    m_FaceVertexNormals.push_back(n[k][2]);
                }
            }

            if(attrib.texcoords.size() > 0) {
                Vec3f tex[3];
                int   t0 = idx0.texcoord_index;
                int   t1 = idx1.texcoord_index;
                int   t2 = idx2.texcoord_index;
                assert(t0 >= 0);
                assert(t1 >= 0);
                assert(t2 >= 0);

                for(int k = 0; k < 2; ++k) {
                    tex[0][k] = attrib.texcoords[2 * t0 + k];
                    tex[1][k] = attrib.texcoords[2 * t1 + k];
                    tex[2][k] = attrib.texcoords[2 * t1 + k];
                }

                for(int k = 0; k < 3; ++k) {
                    m_FaceVertexTexCoord2D.push_back(tex[k][0]);
                    m_FaceVertexTexCoord2D.push_back(tex[k][1]);
                }
            }
        }         // end process current shape
    }

    ////////////////////////////////////////////////////////////////////////////////
    return result;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void MeshLoader::computeFaceVertexData()
{
    if(m_FaceVertexNormals.size() != m_FaceVertices.size()) {
        m_FaceVertexNormals.assign(m_FaceVertices.size(), 0);
        m_FaceVertexColors.assign(m_FaceVertices.size(), 0);
        m_Normals.assign(m_Vertices.size(), 0);

        for(size_t f = 0, f_end = getNFaces(); f < f_end; ++f) {
            // Get index of vertices for the current face
            UInt v0 = m_Faces[3 * f];
            UInt v1 = m_Faces[3 * f + 1];
            UInt v2 = m_Faces[3 * f + 2];

            Vec3f v[3];
            for(Int k = 0; k < 3; ++k) {
                v[0][k] = m_Vertices[3 * v0 + k];
                v[1][k] = m_Vertices[3 * v1 + k];
                v[2][k] = m_Vertices[3 * v2 + k];
            }

            Vec3f faceNormal = glm::normalize(glm::cross(v[1] - v[0], v[2] - v[0]));

            for(Int k = 0; k < 3; ++k) {
                m_Normals[v0 * 3 + k] += faceNormal[k];
                m_Normals[v1 * 3 + k] += faceNormal[k];
                m_Normals[v2 * 3 + k] += faceNormal[k];
            }
        }

        for(size_t f = 0, f_end = getNFaces(); f < f_end; ++f) {
            UInt v0 = m_Faces[3 * f];
            UInt v1 = m_Faces[3 * f + 1];
            UInt v2 = m_Faces[3 * f + 2];

            StdVT_Vec3f fNormals(3, Vec3f(0));
            for(Int k = 0; k < 3; ++k) {
                fNormals[0][k] = m_Normals[3 * v0 + k];
                fNormals[1][k] = m_Normals[3 * v1 + k];
                fNormals[2][k] = m_Normals[3 * v2 + k];
            }

            for(Int k = 0; k < 3; ++k) {
                for(Int l = 0; l < 3; ++l) {
                    m_FaceVertexNormals[9 * f + 3 * k + l] = fNormals[k][l];
                    m_FaceVertexColors[9 * f + 3 * k + l]  = fNormals[k][l];
                }
            }
        }
    }
}
