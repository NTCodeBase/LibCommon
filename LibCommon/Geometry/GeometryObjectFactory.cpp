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

#include <LibCommon/Geometry/GeometryObjectFactory.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NTCodeBase {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
SharedPtr<GeometryObject<N, Real_t>> createGeometry(const String& geometryType, const JParams& jParams) {
    if(geometryType == "Box" || geometryType == "box" || geometryType == "BOX") {
        return std::make_shared<BoxObject<N, Real_t>>(jParams);
    }

    if(geometryType == "Sphere" || geometryType == "sphere" || geometryType == "SPHERE") {
        return std::make_shared<SphereObject<N, Real_t>>(jParams);
    }

    if(geometryType == "Torus" || geometryType == "torus" || geometryType == "TORUS") {
        return std::make_shared<TorusObject<N, Real_t>>(jParams);
    }

    if(geometryType == "Torus28" || geometryType == "torus28" || geometryType == "TORUS28") {
        return std::make_shared<Torus28Object<N, Real_t>>(jParams);
    }

    if(geometryType == "Torus2Inf" || geometryType == "torus2inf" || geometryType == "TORUS2INF") {
        return std::make_shared<Torus2InfObject<N, Real_t>>(jParams);
    }

    if(geometryType == "Torus88" || geometryType == "torus88" || geometryType == "TORUS88") {
        return std::make_shared<Torus88Object<N, Real_t>>(jParams);
    }

    if(geometryType == "TorusInfInf" || geometryType == "torusinfinf" || geometryType == "TORUSINFINF") {
        return std::make_shared<TorusInfInfObject<N, Real_t>>(jParams);
    }

    if(geometryType == "Cylinder" || geometryType == "cylinder" || geometryType == "CYLINDER") {
        return std::make_shared<CylinderObject<N, Real_t>>(jParams);
    }

    if(geometryType == "Cone" || geometryType == "cone" || geometryType == "CONE") {
        return std::make_shared<ConeObject<N, Real_t>>(jParams);
    }

    if(geometryType == "Plane" || geometryType == "plane" || geometryType == "PLANE") {
        return std::make_shared<PlaneObject<N, Real_t>>(jParams);
    }

    if(geometryType == "Triangle" || geometryType == "triangle" || geometryType == "TRIANGLE") {
        return std::make_shared<TriangleObject<N, Real_t>>(jParams);
    }

    if(geometryType == "Hexagon" || geometryType == "hexagon" || geometryType == "HEXAGON") {
        return std::make_shared<HexagonObject<N, Real_t>>(jParams);
    }

    if(geometryType == "TriangularPrism" || geometryType == "triangularprism" || geometryType == "TRIANGULARPRISM") {
        return std::make_shared<TriangularPrismObject<N, Real_t>>(jParams);
    }

    if(geometryType == "HexagonalPrism" || geometryType == "hexagonalprism" || geometryType == "HEXAGONALPRISM") {
        return std::make_shared<HexagonalPrismObject<N, Real_t>>(jParams);
    }

    if(geometryType == "Capsule" || geometryType == "capsule" || geometryType == "CAPSULE") {
        return std::make_shared<CapsuleObject<N, Real_t>>(jParams);
    }

    if(geometryType == "Ellipsoid" || geometryType == "ellipsoid" || geometryType == "ELLIPSOID") {
        return std::make_shared<EllipsoidObject<N, Real_t>>(jParams);
    }

    if(geometryType == "Mesh" || geometryType == "mesh" || geometryType == "MESH" ||
       geometryType == "TriMesh" || geometryType == "trimesh" || geometryType == "TRIMESH") {
        return std::make_shared<TriMeshObject<N, Real_t>>(jParams);
    }

    if(geometryType == "CSGObject" || geometryType == "csgobject" || geometryType == "CSGOBJECT") {
        return std::make_shared<CSGObject<N, Real_t>>(jParams);
    }

    return nullptr;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
SharedPtr<GeometryObject<N, Real_t>> combineGeometryObjects(const StdVT<SharedPtr<GeometryObject<N, Real_t>>>& geometryObjs) {
    if(geometryObjs.size() == 0) {
        return nullptr;
    }
    if(geometryObjs.size() == 1) {
        return geometryObjs.front();
    }

    SharedPtr<CSGObject<N, Real_t>> csgObj = std::make_shared<CSGObject<N, Real_t>>(JParams());

    for(auto& obj : geometryObjs) {
        csgObj->addObject(obj);
    }

    return csgObj;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template SharedPtr<GeometryObject<2, float>> createGeometry(const String& geometryType, const JParams& jParams);
template SharedPtr<GeometryObject<3, float>> createGeometry(const String& geometryType, const JParams& jParams);

template SharedPtr<GeometryObject<2, float>> combineGeometryObjects(const StdVT<SharedPtr<GeometryObject<2, float>>>& geometryObjs);
template SharedPtr<GeometryObject<3, float>> combineGeometryObjects(const StdVT<SharedPtr<GeometryObject<3, float>>>& geometryObjs);

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#ifdef __NT_SUPPORT_DOUBLE_NUMBER
template SharedPtr<GeometryObject<2, double>> createGeometry(const String& geometryType, const JParams& jParams);
template SharedPtr<GeometryObject<3, double>> createGeometry(const String& geometryType, const JParams& jParams);

template SharedPtr<GeometryObject<2, double>> combineGeometryObjects(const StdVT<SharedPtr<GeometryObject<2, double>>>& geometryObjs);
template SharedPtr<GeometryObject<3, double>> combineGeometryObjects(const StdVT<SharedPtr<GeometryObject<3, double>>>& geometryObjs);
#endif
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace NTCodeBase
