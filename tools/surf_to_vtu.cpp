#include <string>
#include <sstream>

#include <tclap/CmdLine.h>
#include "pugixml.hpp"

#include "viennagrid/viennagrid.hpp"
#include "viennagrid/algorithm/inclusion.hpp"
#include "viennagrid/io/vtk_writer.hpp"

#include "viennameshpp/algorithm_pipeline.hpp"

const double eps = 1e-6;


typedef viennagrid::mesh                                          MeshType;
typedef viennagrid::result_of::region<MeshType>::type             RegionType;

typedef viennagrid::result_of::element<MeshType>::type            ElementType;
typedef viennagrid::result_of::point<MeshType>::type              PointType;

void add_svg_polyline(MeshType & mesh, std::string const & polyline, PointType translate)
{
  std::stringstream ss;
  ss << polyline;

  std::string tmp;
  ElementType first_vtx;
  ElementType prev_vtx;
  PointType prev_pt;
  char last_command;


  while (!ss.eof())
  {
    ss >> tmp;

    if (tmp == "m" || tmp == "M" || tmp == "L" || tmp == "l")
    {
      last_command = tmp[0];
      continue;
    }

    ElementType vtx;
    PointType pt;

    if (tmp == "z")
    {
      vtx = first_vtx;
      pt = viennagrid::get_point(first_vtx);
    }
    else
    {
      double x = atof(tmp.substr(0, tmp.find(',')).c_str());
      double y = atof(tmp.substr(tmp.find(',')+1).c_str());

      pt = viennagrid::make_point(x, y);
      if ((last_command == 'm' || last_command == 'l') && !prev_pt.empty())
      {
        pt += prev_pt;
      }
      else
      {
        pt += translate;
      }

      vtx = viennagrid::make_unique_vertex(mesh, pt, eps);

      if (!first_vtx.valid())
      {
        first_vtx = vtx;
      }
    }

    if (prev_vtx.valid())
    {
      viennagrid::make_line(mesh, prev_vtx, vtx);
    }

    prev_vtx = vtx;
    prev_pt = pt;
  }
}

int eliminate_nonconformities(MeshType const & input, MeshType & output)
{
  int counter = 0;

  viennagrid::element_copy_map<double> cpy(output, eps, false);

  typedef viennagrid::result_of::vertex_range<MeshType>::type     VertexRange;
  typedef viennagrid::result_of::iterator<VertexRange>::type        VertexIterator;

  typedef viennagrid::result_of::element_range<MeshType, 1>::type       EdgeRange;
  typedef viennagrid::result_of::iterator<EdgeRange>::type          EdgeIterator;

  EdgeRange edges(input);
  VertexRange vertices(input);

  for (EdgeIterator eit = edges.begin(); eit != edges.end(); ++eit)
  {
    bool found = false;

    for (VertexIterator vit = vertices.begin(); vit != vertices.end(); ++vit)
    {
      if ( (*vit == viennagrid::vertices(*eit)[0]) || (*vit == viennagrid::vertices(*eit)[1]) )
        continue;

      if (viennagrid::is_inside(*eit, viennagrid::get_point(*vit), eps))
      {
        ++counter;
        found = true;

        ElementType ev0 = cpy(viennagrid::vertices(*eit)[0]);
        ElementType ev1 = cpy(viennagrid::vertices(*eit)[1]);
        ElementType v = cpy(*vit);

        try
        {
          viennagrid::make_line(output, ev0, v);
        }
        catch (...) {}

        try
        {
          viennagrid::make_line(output, v, ev1);
        }
        catch (...) {}
      }
    }

    if (!found)
    {
      try
      {
        cpy(*eit);
      }
      catch (...) {}
    }
  }

  return counter;
}


int main(int argc, char **argv)
{
  pugi::xml_document svg_xml;
  TCLAP::UnlabeledValueArg<std::string> surf_filename( "surf_filename", "SURF file name", true, "", "PipelineFile"  );
  TCLAP::UnlabeledValueArg<std::string> vtu_filename( "vtu_filename", "VTU file name", true, "", "PipelineFile"  );

  try
  {
    TCLAP::CmdLine cmd("surf to vtu converter", ' ', "1.0");

    cmd.add( surf_filename );
    cmd.add( vtu_filename );
    cmd.parse( argc, argv );
  }
  catch (TCLAP::ArgException &e)  // catch any exceptions
  {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
  }


  std::ifstream file( surf_filename.getValue().c_str() );
  std::string tmp0, tmp1, tmp2;

  file >> tmp0;
  if (tmp0 != "surfacemesh")
    return -1;

  file >> tmp0;
  int vertex_count = boost::lexical_cast<int>(tmp0);


  MeshType mesh;

  std::vector<ElementType> vertices(vertex_count);

  for (int i = 0; i != vertex_count; ++i)
  {
    file >> tmp0 >> tmp1 >> tmp2;
    PointType pt = viennagrid::make_point( boost::lexical_cast<double>(tmp0),
                                           boost::lexical_cast<double>(tmp1),
                                           boost::lexical_cast<double>(tmp2) );

    vertices[i] = viennagrid::make_vertex(mesh, pt);
  }

  file >> tmp0;
  int triangle_count = boost::lexical_cast<int>(tmp0);

  for (int i = 0; i != triangle_count; ++i)
  {
    file >> tmp0 >> tmp1 >> tmp2;
    ElementType v0 = vertices[boost::lexical_cast<int>(tmp0)-1];
    ElementType v1 = vertices[boost::lexical_cast<int>(tmp1)-1];
    ElementType v2 = vertices[boost::lexical_cast<int>(tmp2)-1];

    viennagrid::make_triangle(mesh, v0, v1, v2);
  }

  std::string output_filename = vtu_filename.getValue().c_str();
  if (output_filename.find(".vtu"))
    output_filename = output_filename.substr( 0, output_filename.find(".vtu") );
  if (output_filename.find(".pvd"))
    output_filename = output_filename.substr( 0, output_filename.find(".pvd") );

  viennagrid::io::vtk_writer<MeshType> writer;
  writer(mesh, output_filename);
}
