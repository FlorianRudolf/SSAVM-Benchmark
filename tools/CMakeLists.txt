cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
include ($ENV{VIENNAUTILSPATH}/cmake/vienna_build_common.cmake NO_POLICY_SCOPE)

vienna_build_setup_project(diss viennamesh viennautils viennagrid)

option(ENABLE_PLUGIN_NETGEN "Enable netgen plugin" OFF)
option(ENABLE_PLUGIN_CGAL "Enable cgal plugin" OFF)
option(ENABLE_PLUGIN_POISSON "Enable poisson plugin" OFF)
option(ENABLE_PLUGIN_METIS "Enable metis plugin" OFF)
option(ENABLE_PLUGIN_TDR "Enable TDR plugin" OFF)
option(ENABLE_PLUGIN_SYMMETRY "Enable symmetry plugin" OFF)
option(ENABLE_PLUGIN_STATISTICS "Enable symmetry plugin" OFF)
option(ENABLE_PLUGIN_MESH_HEALING "Enable mesh healing plugin" OFF)
OPTION(BUILD_SHARED_LIBS  "Build all libraries as shared (dynamic) libraries" ON )
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -felide-constructors -funroll-loops -fstrict-aliasing -Wdisabled-optimization")

vienna_build_add_subprojects()

include_directories($ENV{VIENNAMESHPATH}/external/pugixml-1.5/src)
include_directories($ENV{VIENNAMESHPATH}/external/tclap-1.2.1/include)
# message($ENV{VIENNAMESHPATH}/external)

add_definitions(-DREAL=double -DVOID=int -DNO_TIMER -DTRILIBRARY -DANSI_DECLARATORS)
add_library(diss_common STATIC common.cpp templated_mesh.cpp benchmark.cpp triangle.cpp 3d_geometry.cpp matrix.cpp)
target_link_libraries(diss_common viennameshpp)

add_executable(svg_to_vtu svg_to_vtu.cpp)
target_link_libraries(svg_to_vtu viennameshpp)

add_executable(surf_to_vtu surf_to_vtu.cpp)
target_link_libraries(surf_to_vtu viennameshpp)

add_executable(make_circle_meshes make_circle_meshes.cpp)
target_link_libraries(make_circle_meshes viennameshpp)

add_executable(bench_reflection bench_reflection.cpp)
target_link_libraries(bench_reflection viennameshpp diss_common)

add_executable(bench_rotational_2d_circle bench_rotational_2d_circle.cpp)
target_link_libraries(bench_rotational_2d_circle viennameshpp diss_common)


add_executable(bench_tsv bench_tsv.cpp)
target_link_libraries(bench_tsv viennameshpp diss_common)

add_executable(bench_multi_tsv bench_multi_tsv.cpp)
target_link_libraries(bench_multi_tsv viennameshpp diss_common)

add_executable(bench_bridge bench_bridge.cpp)
target_link_libraries(bench_bridge viennameshpp diss_common)

add_executable(gummel_star gummel_star.cpp)
target_link_libraries(gummel_star viennameshpp diss_common)

add_executable(quantity_difference quantity_difference.cpp)
target_link_libraries(quantity_difference viennameshpp diss_common)


add_executable(compare_solutions compare_solutions.cpp)
target_link_libraries(compare_solutions viennameshpp diss_common)

