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

#include <LibCommon/Geometry/GeometryObjects.h>
#include <LibCommon/Geometry/GeometryObjectFactory.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NTCodeBase {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
bool GeometryObjectFactory<N, Real_t>::registerGeometry(const String& geometryName, CreationFuncPtr creationFunc) {
#ifdef __NT_DEBUG__
    printf("[%s]: Register: %s\n", factoryName().c_str(), geometryName.c_str());
    fflush(stdout);
#endif
    auto[it, bSuccess] = getCreationFuncPtrs().emplace(geometryName, creationFunc);
    __NT_UNUSED(it);
    return bSuccess;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void GeometryObjectFactory<N, Real_t>::registerGeometries() {
    __NT_REQUIRE(registerGeometry(BoxObject<N, Real_t>::name(), BoxObject<N, Real_t>::createGeometry));
    __NT_REQUIRE(registerGeometry(SphereObject<N, Real_t>::name(), SphereObject<N, Real_t>::createGeometry));
    __NT_REQUIRE(registerGeometry(TorusObject<N, Real_t>::name(), TorusObject<N, Real_t>::createGeometry));
    __NT_REQUIRE(registerGeometry(Torus28Object<N, Real_t>::name(), Torus28Object<N, Real_t>::createGeometry));
    __NT_REQUIRE(registerGeometry(Torus2InfObject<N, Real_t>::name(), Torus2InfObject<N, Real_t>::createGeometry));
    __NT_REQUIRE(registerGeometry(Torus88Object<N, Real_t>::name(), Torus88Object<N, Real_t>::createGeometry));
    __NT_REQUIRE(registerGeometry(TorusInfInfObject<N, Real_t>::name(), TorusInfInfObject<N, Real_t>::createGeometry));
    __NT_REQUIRE(registerGeometry(CylinderObject<N, Real_t>::name(), CylinderObject<N, Real_t>::createGeometry));
    __NT_REQUIRE(registerGeometry(ConeObject<N, Real_t>::name(), ConeObject<N, Real_t>::createGeometry));
    __NT_REQUIRE(registerGeometry(CapsuleObject<N, Real_t>::name(), CapsuleObject<N, Real_t>::createGeometry));
    __NT_REQUIRE(registerGeometry(EllipsoidObject<N, Real_t>::name(), EllipsoidObject<N, Real_t>::createGeometry));
    __NT_REQUIRE(registerGeometry(CSGObject<N, Real_t>::name(), CSGObject<N, Real_t>::createGeometry));

    if constexpr(N == 2) {
        __NT_REQUIRE(registerGeometry(PlaneObject<N, Real_t>::name(), PlaneObject<N, Real_t>::createGeometry));
        __NT_REQUIRE(registerGeometry(TriangleObject<N, Real_t>::name(), TriangleObject<N, Real_t>::createGeometry));
        __NT_REQUIRE(registerGeometry(HexagonObject<N, Real_t>::name(), HexagonObject<N, Real_t>::createGeometry));
    } else {
        __NT_REQUIRE(registerGeometry(TriangularPrismObject<N, Real_t>::name(), TriangularPrismObject<N, Real_t>::createGeometry));
        __NT_REQUIRE(registerGeometry(HexagonalPrismObject<N, Real_t>::name(), HexagonalPrismObject<N, Real_t>::createGeometry));
        __NT_REQUIRE(registerGeometry(TriMeshObject<N, Real_t>::name(), TriMeshObject<N, Real_t>::createGeometry));
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
SharedPtr<GeometryObject<N, Real_t>> GeometryObjectFactory<N, Real_t>::createGeometry(const String& geometryName, const JParams& jParams) {
    if(!s_bInitialized) {
        registerGeometries();
        s_bInitialized = true;
    }
    ////////////////////////////////////////////////////////////////////////////////
    auto& creationFuncs = getCreationFuncPtrs();
    __NT_REQUIRE(creationFuncs.find(geometryName) != creationFuncs.end());
    auto geoObj = creationFuncs.at(geometryName)(jParams); // call the creationFunc
    __NT_REQUIRE(geoObj != nullptr);
    return geoObj;
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
template<Int N, class Real_t>
std::unordered_map<String, typename GeometryObjectFactory<N, Real_t>::CreationFuncPtr>&
GeometryObjectFactory<N, Real_t>::getCreationFuncPtrs() {
    static std::unordered_map<String, CreationFuncPtr> creationFuncPtrs;
    return creationFuncPtrs;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
String GeometryObjectFactory<N, Real_t>::factoryName() {
    static String name = String("GeometryObjectFactory-") + std::to_string(N) + String("D-") +
                         (std::is_same_v<Real_t, float> ? String("float") : String("double"));
    return name;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
__NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(GeometryObjectFactory)
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace NTCodeBase
