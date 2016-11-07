#ifndef _VTK_READER_H_
#define _VTK_READER_H_

#include "pugixml-1.5/src/pugixml.hpp"



/* creates a vertex in a mesh, the ID of the created vertex will be returned in vertex_id (optional, will be ignored if vertex_id is NULL). mesh has to be a root mesh which doesn't have any child meshes. */
VIENNAGRID_DYNAMIC_EXPORT viennagrid_error viennagrid_mesh_vertex_create_unique(viennagrid_mesh mesh,
                                                                                const viennagrid_numeric * coords,
                                                                                viennagrid_numeric tolerance,
                                                                                viennagrid_element_id * vertex_id)
{
  viennagrid_dimension geometric_dimension;
  viennagrid_mesh_geometric_dimension_get(mesh, &geometric_dimension);

  viennagrid_element_id * vertex_ids_begin;
  viennagrid_element_id * vertex_ids_end;
  viennagrid_mesh_elements_get(mesh, 0, &vertex_ids_begin, &vertex_ids_end);

  for (viennagrid_element_id * vit = vertex_ids_begin; vit != vertex_ids_end; ++vit)
  {
    viennagrid_numeric * tmp;
    viennagrid_mesh_vertex_coords_get(mesh, *vit, &tmp);

    viennagrid_bool result;

    viennagrid_point_is_equal(geometric_dimension, coords, tmp, tolerance, &result);

    if (result == VIENNAGRID_TRUE)
    {
      *vertex_id = *vit;
      return VIENNAGRID_SUCCESS;
    }
  }

  return viennagrid_mesh_vertex_create(mesh, coords, vertex_id);
}


template<typename T>
std::vector<T> string_to_vector(char const * str)
{
  std::stringstream pds(str);
  std::vector<T> results;
  while (!pds.eof())
  {
    viennagrid_numeric tmp;
    pds >> tmp;
    results.push_back(tmp);
  }
  results.pop_back();

  return results;
}


VIENNAGRID_DYNAMIC_EXPORT viennagrid_error viennagrid_mesh_io_add_vtu(viennagrid_mesh_io mesh_io,
                                                                      const char * filename,
                                                                      viennagrid_region region)
{
  viennagrid_numeric tolerance = 1e-8;

  viennagrid_mesh mesh;
  viennagrid_mesh_io_mesh_get(mesh_io, &mesh);
  if (!mesh)
  {
    viennagrid_mesh_create(&mesh);
    viennagrid_mesh_property_set(mesh, VIENNAGRID_PROPERTY_BOUNDARY_LAYOUT, VIENNAGRID_BOUNDARY_LAYOUT_FULL);
    viennagrid_mesh_io_mesh_set(mesh_io, mesh);
  }


  pugi::xml_document VTUFile;
  pugi::xml_parse_result result = VTUFile.load_file(filename);

  pugi::xml_node VTKFile = VTUFile.child("VTKFile");
  pugi::xml_node UnstructuredGrid = VTKFile.child("UnstructuredGrid");
  pugi::xml_node Piece = UnstructuredGrid.child("Piece");


  pugi::xml_attribute num_points_attribute = Piece.attribute("NumberOfPoints");
  pugi::xml_attribute num_cells_attribute = Piece.attribute("NumberOfCells");

  int num_points = num_points_attribute.as_int();
  int num_cells = num_cells_attribute.as_int();


  std::vector<viennagrid_element_id> vertex_mapping;

  pugi::xml_node Points = Piece.child("Points");
  {
    pugi::xml_node DataArray = Points.child("DataArray");

    pugi::xml_attribute type_attribute = DataArray.attribute("type");
    pugi::xml_attribute NumberOfComponents_attribute = DataArray.attribute("NumberOfComponents");
    pugi::xml_attribute format_attribute = DataArray.attribute("format");

//     if (std::string(type_attribute.as_string()) != "Float64")
//       return VIENNAGRID_ERROR_UNSPECIFIED_ERROR;
//     if (std::string(format_attribute.as_string()) != "ascii")
//       return VIENNAGRID_ERROR_UNSPECIFIED_ERROR;

    int file_geometric_dimension = NumberOfComponents_attribute.as_int();
    int geometric_dimension = 0;

    std::vector<viennagrid_numeric> coords = string_to_vector<viennagrid_numeric>(DataArray.child_value());
    for (int i = 0; i != num_points; ++i)
    {
      viennagrid_numeric * point = &coords[i*file_geometric_dimension];
      int tmp = file_geometric_dimension-1;
      for (; tmp >= 0; --tmp)
      {
        if (std::abs(point[tmp]) >= tolerance)
          break;
      }

//       std::cout << "TMP " << tmp << std::endl;
      geometric_dimension = std::max(geometric_dimension, tmp);
    }

    ++geometric_dimension;
//     std::cout << "GEO DIM " << geometric_dimension << std::endl;

    {
      viennagrid_dimension mesh_geometric_dimension;
      viennagrid_mesh_geometric_dimension_get(mesh, &mesh_geometric_dimension);
      if (mesh_geometric_dimension <= 0) // if mesh is not initialized, set the geometric dimension
        viennagrid_mesh_geometric_dimension_set(mesh, geometric_dimension);
      else if (mesh_geometric_dimension != geometric_dimension)
        return VIENNAGRID_ERROR_UNSPECIFIED_ERROR;
    }



    vertex_mapping.resize(num_points);
    for (int i = 0; i != num_points; ++i)
      viennagrid_mesh_vertex_create_unique(mesh, &coords[i*file_geometric_dimension], tolerance, &vertex_mapping[i]);
  }


  pugi::xml_node Cells = Piece.child("Cells");
  {
    pugi::xml_node ConnectivityDataArray;
    pugi::xml_node OffsetsDataArray;
    pugi::xml_node TypesDataArray;

    {
      pugi::xml_node DataArray = Cells.child("DataArray");
      do
      {
        pugi::xml_attribute Name_attribute = DataArray.attribute("Name");

        if (std::string(Name_attribute.as_string()) == "connectivity")
          ConnectivityDataArray = DataArray;
        else if (std::string(Name_attribute.as_string()) == "offsets")
          OffsetsDataArray = DataArray;
        else if (std::string(Name_attribute.as_string()) == "types")
          TypesDataArray = DataArray;

        DataArray = DataArray.next_sibling("DataArray");
      }
      while (DataArray);

      std::vector<viennagrid_int> connectivity = string_to_vector<viennagrid_int>(ConnectivityDataArray.child_value());
      std::vector<viennagrid_int> offsets = string_to_vector<viennagrid_int>(OffsetsDataArray.child_value());
      std::vector<viennagrid_element_type> types;

      {
        std::vector<viennagrid_int> tmp = string_to_vector<viennagrid_int>(TypesDataArray.child_value());
        types.resize(tmp.size());
        for (std::size_t i = 0; i != tmp.size(); ++i)
        {
          switch (tmp[i])
          {
            case 1:
              types[i] = VIENNAGRID_ELEMENT_TYPE_VERTEX;
              break;
            case 3:
              types[i] = VIENNAGRID_ELEMENT_TYPE_LINE;
              break;
            case 5:
              types[i] = VIENNAGRID_ELEMENT_TYPE_TRIANGLE;
              break;
            case 7:
              types[i] = VIENNAGRID_ELEMENT_TYPE_POLYGON;
              break;
            case 9:
              types[i] = VIENNAGRID_ELEMENT_TYPE_QUADRILATERAL;
              break;
            case 10:
              types[i] = VIENNAGRID_ELEMENT_TYPE_TETRAHEDRON;
              break;
            case 12:
              types[i] = VIENNAGRID_ELEMENT_TYPE_HEXAHEDRON;
              break;

            default:
              return VIENNAGRID_ERROR_UNSPECIFIED_ERROR;
          }
        }
      }

      std::vector<viennagrid_region_id> regions;
      if (region)
      {
        viennagrid_region_id region_id;
        viennagrid_region_id_get(region, &region_id);

        regions.resize(types.size(), region_id);
      }

      for (std::size_t i = 0; i != connectivity.size(); ++i)
      {
        connectivity[i] = vertex_mapping[ connectivity[i] ];
      }

      offsets.insert(offsets.begin(), 0);

      viennagrid_element_id first_cell;
      viennagrid_mesh_element_batch_create(mesh, types.size(), &types[0], &offsets[0], &connectivity[0], 0, &first_cell);
    }
  }

  viennagrid_dimension cell_dimension;
  viennagrid_mesh_cell_dimension_get(mesh, &cell_dimension);


  pugi::xml_node PointData = Piece.child("PointData");
  if (PointData)
  {
    pugi::xml_node DataArray = PointData.child("DataArray");
    do
    {
      pugi::xml_attribute Name_attribute = DataArray.attribute("Name");
      std::string data_name = Name_attribute.as_string();

      std::vector<viennagrid_numeric> data = string_to_vector<viennagrid_numeric>(DataArray.child_value());

      viennagrid_quantity_field qf;
      viennagrid_mesh_io_quantity_field_get(mesh_io, data_name.c_str(), &qf);
      if (!qf)
      {
        viennagrid_quantity_field_create(&qf);
        viennagrid_quantity_field_name_set(qf, data_name.c_str());
        viennagrid_quantity_field_init(qf, 0, VIENNAGRID_QUANTITY_FIELD_TYPE_NUMERIC, 1, VIENNAGRID_QUANTITY_FIELD_STORAGE_DENSE);
        viennagrid_mesh_io_quantity_field_set(mesh_io, qf);
      }

      for (int i = 0; i != num_points; ++i)
        viennagrid_quantity_field_value_set(qf, vertex_mapping[i], &data[i]);

      DataArray = DataArray.next_sibling("DataArray");
    }
    while (DataArray);
  }



  pugi::xml_node CellData = Piece.child("CellData");
  if (CellData)
  {
    pugi::xml_node DataArray = CellData.child("DataArray");
    while (DataArray)
    {
      pugi::xml_attribute Name_attribute = DataArray.attribute("Name");
      std::string data_name = Name_attribute.as_string();

      std::vector<viennagrid_numeric> data = string_to_vector<viennagrid_numeric>(DataArray.child_value());

      viennagrid_quantity_field qf;
      viennagrid_mesh_io_quantity_field_get(mesh_io, data_name.c_str(), &qf);
      if (!qf)
      {
        viennagrid_quantity_field_create(&qf);
        viennagrid_quantity_field_name_set(qf, data_name.c_str());
        viennagrid_quantity_field_init(qf, cell_dimension, VIENNAGRID_QUANTITY_FIELD_TYPE_NUMERIC, 1, VIENNAGRID_QUANTITY_FIELD_STORAGE_DENSE);
        viennagrid_mesh_io_quantity_field_set(mesh_io, qf);
      }

      for (int i = 0; i != num_points; ++i)
        viennagrid_quantity_field_value_set(qf, vertex_mapping[i], &data[i]);

      DataArray = DataArray.next_sibling("DataArray");
    }
  }

  return VIENNAGRID_SUCCESS;
}



VIENNAGRID_DYNAMIC_EXPORT viennagrid_error viennagrid_mesh_io_read_pvd(viennagrid_mesh_io mesh_io,
                                                                       const char * fn)
{
  std::string filename = fn;
  std::string path = filename.substr(0, filename.rfind("/"));

  viennagrid_mesh mesh;
  viennagrid_mesh_io_mesh_get(mesh_io, &mesh);
  if (!mesh)
  {
    viennagrid_mesh_create(&mesh);
    viennagrid_mesh_property_set(mesh, VIENNAGRID_PROPERTY_BOUNDARY_LAYOUT, VIENNAGRID_BOUNDARY_LAYOUT_FULL);
    viennagrid_mesh_io_mesh_set(mesh_io, mesh);
  }

  pugi::xml_document PVDFile;
  pugi::xml_parse_result result = PVDFile.load_file(filename.c_str());

  pugi::xml_node VTKFile = PVDFile.child("VTKFile");
  pugi::xml_attribute type_attribute = VTKFile.attribute("type");

  if (std::string(type_attribute.as_string()) != "Collection")
    return VIENNAGRID_ERROR_UNSPECIFIED_ERROR;

  pugi::xml_node Collection = VTKFile.child("Collection");

  pugi::xml_node DataSet = Collection.child("DataSet");
  while (DataSet)
  {
    pugi::xml_attribute part_attribute = DataSet.attribute("part");
    pugi::xml_attribute file_attribute = DataSet.attribute("file");

    viennagrid_region_id region_id = part_attribute.as_int();
    std::string vtu_filename = path + "/" + file_attribute.as_string();

    viennagrid_region region;
    viennagrid_mesh_region_get_or_create(mesh, region_id, &region);

//     std::cout << "Add VTU " << vtu_filename << std::endl;
    viennagrid_mesh_io_add_vtu(mesh_io, vtu_filename.c_str(), region);

    DataSet = DataSet.next_sibling("DataSet");
  }

  return VIENNAGRID_SUCCESS;
}


#endif
