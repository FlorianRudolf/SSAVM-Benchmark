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


PointType rotate(PointType const & pt, double angle)
{
  double ca = std::cos(angle);
  double sa = std::sin(angle);

  return viennagrid::make_point( ca*pt[0] - sa*pt[1], sa*pt[0] + ca*pt[1] );
}


int main(int argc, char **argv)
{
  pugi::xml_document svg_xml;

  TCLAP::ValueArg<int> rotational_symmetry_order("s","rotational-symmetry-order", "Rotational Symmetry Order", true, -1, "int");
  TCLAP::ValueArg<double> radius("r","radius", "Radius", false, 1.0, "double");
  TCLAP::UnlabeledValueArg<std::string> filename( "filename", "File name", true, "", "FileName"  );

  try
  {
    TCLAP::CmdLine cmd("Make circle mesh", ' ', "1.0");

    cmd.add( rotational_symmetry_order );
    cmd.add( radius );
    cmd.add( filename );

    cmd.parse( argc, argv );
  }
  catch (TCLAP::ArgException &e)  // catch any exceptions
  {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
  }


  double angle = 2*M_PI/rotational_symmetry_order.getValue();

  std::string output_filename = filename.getValue().c_str();
  if (output_filename.find(".vtu"))
    output_filename = output_filename.substr( 0, output_filename.find(".vtu") );
  if (output_filename.find(".pvd"))
    output_filename = output_filename.substr( 0, output_filename.find(".pvd") );


    {
      MeshType mesh;

      PointType first_pt = viennagrid::make_point(0, radius.getValue());
      ElementType first_vtx = viennagrid::make_vertex(mesh, first_pt);

      ElementType prev_vtx = first_vtx;

      for (int i = 1; i != rotational_symmetry_order.getValue(); ++i)
      {
        PointType pt = rotate( first_pt, i*angle );
        ElementType vtx = viennagrid::make_vertex(mesh, pt);

        viennagrid::make_line(mesh, prev_vtx, vtx);

        prev_vtx = vtx;
      }

      viennagrid::make_line(mesh, prev_vtx, first_vtx);



      viennagrid::io::vtk_writer<MeshType> writer;
      writer(mesh, output_filename);
    }



    {
      MeshType mesh;

      PointType top_pt = viennagrid::make_point(0, radius.getValue());
      ElementType top_vtx = viennagrid::make_vertex(mesh, top_pt);

      PointType left_pt = (top_pt + rotate(top_pt, angle)) / 2;
      ElementType left_vtx = viennagrid::make_vertex(mesh, left_pt);

      PointType right_pt = (top_pt + rotate(top_pt, -angle)) / 2;
      ElementType right_vtx = viennagrid::make_vertex(mesh, right_pt);

      PointType center_pt = viennagrid::make_point(0, 0);
      ElementType center_vtx = viennagrid::make_vertex(mesh, center_pt);

      viennagrid::make_line(mesh, top_vtx, left_vtx);
      viennagrid::make_line(mesh, left_vtx, center_vtx);
      viennagrid::make_line(mesh, center_vtx, right_vtx);
      viennagrid::make_line(mesh, right_vtx, top_vtx);

      viennagrid::io::vtk_writer<MeshType> writer;
      writer(mesh, output_filename + "_slice");
    }

    {
      MeshType mesh;

      PointType top_pt = viennagrid::make_point(0, radius.getValue());
      ElementType top_vtx = viennagrid::make_vertex(mesh, top_pt);

      PointType left_pt = (top_pt + rotate(top_pt, angle)) / 2;
      ElementType left_vtx = viennagrid::make_vertex(mesh, left_pt);

      PointType center_pt = viennagrid::make_point(0, 0);
      ElementType center_vtx = viennagrid::make_vertex(mesh, center_pt);

      viennagrid::make_line(mesh, top_vtx, left_vtx);
      viennagrid::make_line(mesh, left_vtx, center_vtx);
      viennagrid::make_line(mesh, center_vtx, top_vtx);

      viennagrid::io::vtk_writer<MeshType> writer;
      writer(mesh, output_filename + "_halfslice");
    }
}
