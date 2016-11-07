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
  TCLAP::UnlabeledValueArg<std::string> svg_filename( "svg_filename", "SVG file name", true, "", "PipelineFile"  );
  TCLAP::UnlabeledValueArg<std::string> vtu_filename( "vtu_filename", "SVG file name", true, "", "PipelineFile"  );

  try
  {
    TCLAP::CmdLine cmd("SVG to vtu converter", ' ', "1.0");

    cmd.add( svg_filename );
    cmd.add( vtu_filename );
    cmd.parse( argc, argv );

    pugi::xml_parse_result result = svg_xml.load_file( svg_filename.getValue().c_str() );

    if (!result)
    {
      std::cerr << "Error loading or parsing XML file " << svg_filename.getValue().c_str() << std::endl;
      std::cerr << "XML error: " << result.description() << std::endl;
      return 0;
    }
  }
  catch (TCLAP::ArgException &e)  // catch any exceptions
  {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
  }


  MeshType mesh;

  pugi::xml_node svg_node = svg_xml.child("svg");
  pugi::xml_node g_node = svg_node.child("g");


  for (pugi::xml_node path_node = g_node.child("path");
        path_node;
        path_node = path_node.next_sibling("path"))
  {
    std::string tmp = path_node.attribute("d").as_string();
    std::cout << tmp << std::endl;

    PointType translate = viennagrid::make_point(0,0);

    pugi::xml_attribute transform_attribute = path_node.attribute("transform");
    if ( !transform_attribute.empty() )
    {
      std::string transform_str = transform_attribute.as_string();
      if (transform_str.find("translate") != std::string::npos)
      {
//         std::cout << "translate x = " << ) << std::endl;
//         std::cout << "translate y = " << transform_str.substr( transform_str.find(",")+1, transform_str.find(")")-transform_str.find(",")-1 ) << std::endl;
        double x = atof( transform_str.substr( transform_str.find("(")+1, transform_str.find(",")-transform_str.find("(")-1 ).c_str() );
        double y = atof( transform_str.substr( transform_str.find(",")+1, transform_str.find(")")-transform_str.find(",")-1 ).c_str() );

        translate = viennagrid::make_point(x,y);
      }
      else
      {
        std::cerr << transform_str << " is an invalid transform (currently only translate is supported)" << std::endl;
        return 0;
      }
    }

    add_svg_polyline(mesh, tmp, translate);
  }




  int fixes = 0;
  do
  {
    MeshType tmp;
    fixes = eliminate_nonconformities(mesh, tmp);
    mesh = tmp;

    std::cout << "Fix count = " << fixes << std::endl;
  } while (fixes != 0);



  std::string output_filename = vtu_filename.getValue().c_str();
  if (output_filename.find(".vtu"))
    output_filename = output_filename.substr( 0, output_filename.find(".vtu") );
  if (output_filename.find(".pvd"))
    output_filename = output_filename.substr( 0, output_filename.find(".pvd") );

  viennagrid::io::vtk_writer<MeshType> writer;
  writer(mesh, output_filename);



//   viennamesh::context_handle context;
//
//   viennamesh::algorithm_handle mesher = context.make_algorithm("triangle_make_mesh");
//   mesher.set_input( "mesh", mesh.internal() );
//   mesher.set_input("cell_size", 10.0);
//   {
//     viennamesh::LoggingStack s("tetgen_make_mesh");
//     mesher.run();
//   }
//
//   viennamesh::algorithm_handle mesh_writer = context.make_algorithm("mesh_writer");
//   mesh_writer.set_default_source(mesher);
//   mesh_writer.set_input( "filename", "half_aircraft_mesh.vtu" );
//   {
//     viennamesh::LoggingStack s("mesh_writer");
//     mesh_writer.run();
//   }
}
