#include <tclap/CmdLine.h>

#include "common.hpp"
#include "templated_mesh.hpp"
#include "benchmark.hpp"
#include "3d_geometry.hpp"

#include "viennagrid/io/vtk_reader.hpp"
#include "viennagrid/algorithm/quantity_interpolate.hpp"


int main(int argc, char **argv)
{
  TCLAP::ValueArg<std::string> highres_filename("k","highres-filename", "Highres mesh filename", false, "", "string");
  TCLAP::ValueArg<std::string> lowres_filename("l","lowres-filename", "Lowres mesh filename", false, "", "string");

  TCLAP::ValueArg<std::string> highres_qfname("q","highres-quantityname", "Highres quantity name", false, "", "string");
  TCLAP::ValueArg<std::string> lowres_qfname("p","lowres-quantityname", "Lowres quantity name", false, "", "string");

  try
  {
    TCLAP::CmdLine cmd("Quantity difference", ' ', "1.0");

    cmd.add( highres_filename );
    cmd.add( lowres_filename );

    cmd.add( highres_qfname );
    cmd.add( lowres_qfname );

    cmd.parse( argc, argv );
  }
  catch (TCLAP::ArgException &e)  // catch any exceptions
  {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return 0;
  }



  MeshType highres_mesh;
  viennagrid::quantity_field highres_quantity_field;
  {
    viennagrid::io::vtk_reader<viennagrid::mesh> reader;
    reader(highres_mesh, highres_filename.getValue());
    std::vector<viennagrid::quantity_field> quantity_fields = reader.quantity_fields();

    for (std::size_t i = 0; i != quantity_fields.size(); ++i)
    {
      if (quantity_fields[i].get_name() == highres_qfname.getValue())
      {
        highres_quantity_field = quantity_fields[i];
        break;
      }
    }
  }



  MeshType lowres_mesh;
  viennagrid::quantity_field lowres_quantity_field;
  {
    viennagrid::io::vtk_reader<viennagrid::mesh> reader;
    reader(lowres_mesh, highres_filename.getValue());
    std::vector<viennagrid::quantity_field> quantity_fields = reader.quantity_fields();

    for (std::size_t i = 0; i != quantity_fields.size(); ++i)
    {
      if (quantity_fields[i].get_name() == lowres_qfname.getValue())
      {
        lowres_quantity_field = quantity_fields[i];
        break;
      }
    }
  }


  viennagrid::quantity_field lowres_quantity_field_interpolated_to_highres( 0, highres_quantity_field.values_per_quantity(), highres_quantity_field.storage_layout() );
  viennagrid::interpolate_vertex_quantity( lowres_mesh, lowres_quantity_field, highres_mesh, lowres_quantity_field_interpolated_to_highres, 0 );

  ElementRangeType vertices(highres_mesh, 0);
  for (ElementRangeIterator vit = vertices.begin(); vit != vertices.end(); ++vit)
  {
    double lowres_value = lowres_quantity_field_interpolated_to_highres.get(*vit);
    double highres_value = highres_quantity_field.get(*vit);
    std::cout << lowres_value << " " << highres_value << std::endl;
  }

  return 0;
}

