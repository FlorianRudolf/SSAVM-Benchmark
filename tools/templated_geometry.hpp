#ifndef TEMPLATED_GEOMETRY_HPP_
#define TEMPLATED_GEOMETRY_HPP_

#include "common.hpp"
#include "pugixml.hpp"


struct TemplatedGeometry;

struct GeometryTemplate
{
  MeshType          boundary_2d;
  viennagrid_plc    boundary_3d;
};

struct GeometryInstance
{
  GeometryInstance(TemplatedGeometry & templated_geometry_) : templated_geometry(&templated_geometry_) {}

  int geometry_template_index;
  int transformation_index;
  int region_id;

  TemplatedGeometry * templated_geometry;
};

struct TemplatedGeometry
{
  std::vector<GeometryTemplate>   templates;
  std::vector<Transformation>     transformations;
  std::vector<GeometryInstance>   instances;

  PointType min;
  PointType max;
};




#endif
