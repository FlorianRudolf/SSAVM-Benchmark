#ifndef DISS_3D_GEOMETRY_HPP_
#define DISS_3D_GEOMETRY_HPP_


#include "triangle.h"
#include "common.hpp"
#include "benchmark.hpp"

struct TriangleWrapper : public diss::triangulateio
{
  TriangleWrapper()
  {
    pointlist = NULL;
    pointattributelist = NULL;
    pointmarkerlist = NULL;
    numberofpoints = 0;
    numberofpointattributes = 0;

    trianglelist = NULL;
    triangleattributelist = NULL;
    trianglearealist = NULL;
    neighborlist = NULL;
    numberoftriangles = 0;
    numberofcorners = 0;
    numberoftriangleattributes = 0;

    segmentlist = NULL;
    segmentmarkerlist = NULL;
    numberofsegments = 0;

    holelist = NULL;
    numberofholes = 0;

    regionlist = NULL;
    numberofregions = 0;

    edgelist = NULL;
    edgemarkerlist = NULL;
    normlist = NULL;
    numberofedges = 0;
  }

  ~TriangleWrapper()
  {
    if (pointlist) free(pointlist);
    if (pointattributelist) free(pointattributelist);
    if (pointmarkerlist) free(pointmarkerlist);

    if (trianglelist) free(trianglelist);
    if (triangleattributelist) free(triangleattributelist);
    if (trianglearealist) free(trianglearealist);
    if (neighborlist) free(neighborlist);

    if (segmentlist) free(segmentlist);
    if (segmentmarkerlist) free(segmentmarkerlist);

    if (edgelist) free(edgelist);
    if (edgemarkerlist) free(edgemarkerlist);
    if (normlist) free(normlist);

    if (holelist) free(holelist);
    if (regionlist) free(regionlist);
  }


  void init_points(int num_points)
  {
    numberofpoints = num_points;
    if (pointlist) free(pointlist);
    pointlist = (REAL*)malloc(sizeof(REAL) * 2 * numberofpoints);
  }

  void init_segments(int num_segments)
  {
    numberofsegments = num_segments;
    if (segmentlist) free(segmentlist);
    segmentlist = (int*)malloc(sizeof(int) * 2 * numberofsegments);
  }

  void init_hole_points(int num_hole_points)
  {
    numberofholes = num_hole_points;
    if (holelist) free(holelist);
    holelist = (REAL*)malloc(sizeof(REAL) * 2 * numberofholes);
  }


  MeshType to_viennagrid() const
  {
    MeshType result;
    std::vector<ElementType> vertices(numberofpoints);

    for (int i = 0; i < numberofpoints; ++i)
    {
      vertices[i] = viennagrid::make_vertex( result, viennagrid::make_point(pointlist[2*i+0], pointlist[2*i+1]) );
    }

    for (int i = 0; i < numberoftriangles; ++i)
    {
      ElementType cell = viennagrid::make_triangle(
        result,
        vertices[ trianglelist[3*i+0] ],
        vertices[ trianglelist[3*i+1] ],
        vertices[ trianglelist[3*i+2] ]
      );
    }

    return result;
  }
};


struct TriangleMesh
{
  TriangleMesh() : generated_(false) {}

  bool generated() const { return generated_; }
  MeshType to_viennagrid() const { return mesh.to_viennagrid(); }

  void generate(std::string const & option_string)
  {
    if (!generated())
    {
      triangulate( const_cast<char*>(option_string.c_str()), &geometry, &mesh, NULL);
      generated_ = true;
    }
  }

  void reset_mesh()
  {
    generated_ = false;
    mesh = TriangleWrapper();
  }

  TriangleWrapper geometry;
  TriangleWrapper mesh;
  bool generated_;
};



struct TemplatedFacet
{
  TemplatedFacet() : facet_template(new TriangleMesh()) {}

  boost::shared_ptr<TriangleMesh> facet_template;
  Transformation from_2d;
};


struct TemplatedSurface
{
  std::vector<TemplatedFacet> facets;

  void reset_meshes()
  {
    for (std::size_t i = 0; i != facets.size(); ++i)
      facets[i].facet_template->reset_mesh();
  }

  void generate_meshes(std::string const & options)
  {
    for (std::size_t i = 0; i != facets.size(); ++i)
    {
      facets[i].facet_template->generate(options);

      facets[i].facet_template->mesh.numberofholes = 0;
      facets[i].facet_template->mesh.holelist = NULL;
    }
  }

  MeshType compose(PointType const & min, PointType const & max, viennagrid_numeric eps) const
  {
    std::vector<MeshType> vgridmesh;
    MeshType surface_mesh;
    CopyMapType cm(surface_mesh, eps, min, max);

    for (std::size_t i = 0; i != facets.size(); ++i)
    {
      MeshType tmpmesh = facets[i].facet_template->to_viennagrid();

      MeshType transformed_3d = facets[i].from_2d(tmpmesh);
      vgridmesh.push_back(transformed_3d);
      cm(transformed_3d);
    }

    return surface_mesh;
  }
};


Transformation find_transform(TriangleWrapper const & from,
                              TriangleWrapper const & to,
                              viennagrid_numeric eps,
                              bool swap_orientation = false);




viennagrid_plc make_rotational_3d_geometry(TSVconfig const & cfg, viennagrid_numeric eps);
viennagrid_plc make_rotational_3d_geometry_nonsym(TSVconfig const & cfg, viennagrid_numeric eps);
viennagrid_plc make_rotational_3d_geometry_nonsym_middle(TSVconfig const & cfg);
viennagrid_plc make_multi_tsv_geometry(TSVconfig const & tsv_cfg, std::vector<PointType> const & centers, viennagrid_numeric eps);
// viennagrid_plc make_beam(BridgeConfig const & cfg, viennagrid_numeric eps);
// viennagrid_plc make_beam_connector(BridgeConfig const & cfg, viennagrid_numeric eps);

#endif
