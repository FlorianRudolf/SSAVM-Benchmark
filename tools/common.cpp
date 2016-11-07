#include "common.hpp"

std::vector<std::string> split_string( std::string const & str, std::string const & delimiter )
{
  std::vector<std::string> tokens;

  std::string::size_type pos = 0;
  while (pos < str.size())
  {
    std::string::size_type new_pos = str.find(delimiter, pos);
    if (new_pos == std::string::npos)
    {
      tokens.push_back( str.substr(pos, str.size()-pos) );
      return tokens;
    }

    tokens.push_back( str.substr(pos, new_pos-pos) );
    pos = new_pos+delimiter.size();
  }

  if (pos == str.size())
    return tokens;

  std::cout << "something went wrong..." << std::endl;

  return std::vector<std::string>();
}

std::vector<std::string> split_string_brackets( std::string const & str, std::string const & delimiter )
{
  std::vector<std::string> tokens;

  std::vector<int> bracket_count(str.size());

  std::string::size_type pos = 0;
  int bc = 0;
  for (pos = 0; pos != str.size(); ++pos)
  {
    if (str[pos] == '(')
      ++bc;

    bracket_count[pos] = bc;

    if (str[pos] == ')')
      --bc;
  }

  pos = 0;
  std::string::size_type old_pos = 0;
  while (pos < str.size())
  {
    std::string::size_type new_pos = str.find(delimiter, pos);
    if (new_pos == std::string::npos)
    {
      tokens.push_back( str.substr(old_pos, str.size()-old_pos) );
      return tokens;
    }

    if (bracket_count[new_pos] != 0)
    {
      pos = new_pos+delimiter.size();
      continue;
    }


    tokens.push_back( str.substr(old_pos, new_pos-old_pos) );
    pos = new_pos+delimiter.size();
    old_pos = pos;
  }

  if (pos == str.size())
    return tokens;

  std::cout << "something went wrong..." << std::endl;

  return std::vector<std::string>();
}



viennamesh::context_handle context()
{
  static viennamesh::context_handle ctx;
  return ctx;
}



void prod(MatrixType const & matrix, viennagrid_numeric const * pt, viennagrid_numeric * result)
{
  for (int row = 0; row != matrix.row_count; ++row)
  {
    result[row] = 0;
    for (int col = 0; col != matrix.column_count; ++col)
    {
      result[row] += pt[col]*matrix(row,col);
    }
  }
}


PointType prod(MatrixType const & matrix, PointType const & pt)
{
  assert( pt.size() == matrix.column_count );
  PointType result( matrix.row_count );
  prod(matrix, &pt[0], &result[0]);
  return result;
}

MatrixType prod(MatrixType const & A1, MatrixType const & A2)
{
  assert( A1.column_count == A2.row_count );
  MatrixType result(A1.row_count, A2.column_count);

  for (int row = 0; row != result.row_count; ++row)
  {
    for (int column = 0; column != result.column_count; ++column)
    {
      result(row,column) = 0;
      for (int i = 0; i != A1.column_count; ++i)
        result(row,column) += A1(row,i) * A2(i,column);
    }
  }

  return result;
}




// std::pair<Transformation, Transformation> Transformation::make_facet_projection(viennagrid_plc plc, viennagrid_element_id facet_id, viennagrid_numeric eps)
// {
//   viennagrid_dimension geo_dim;
//   viennagrid_plc_geometric_dimension_get(plc, &geo_dim);
//
//   viennagrid_element_id * facet_vertices_begin;
//   viennagrid_element_id * facet_vertices_end;
//   viennagrid_plc_boundary_elements(plc, facet_id, 0, &facet_vertices_begin, &facet_vertices_end);
//   viennagrid_int facet_vertex_count = facet_vertices_end-facet_vertices_begin;
//
//   PointType pt0;
//   PointType pt1;
//   PointType pt2;
//
//   viennagrid_element_id * vit0 = facet_vertices_begin;
//   viennagrid_element_id * vit1 = vit0; ++vit1;
//   viennagrid_element_id * vit2 = vit1; ++vit2;
//
//   bool found;
//   for (; vit0 != facet_vertices_end; ++vit0)
//   {
//     pt0 = get_point(plc, *vit0);
//     for (; vit1 != facet_vertices_end; ++vit1)
//     {
//       pt1 = get_point(plc, *vit1);
//       for (; vit2 != facet_vertices_end; ++vit2)
//       {
//         pt2 = get_point(plc, *vit2);
//
//         PointType c = viennagrid::cross_prod(pt1-pt0, pt2-pt0);
//         if (viennagrid::norm_2(c) > eps)
//         {
//           found = true;
//           break;
//         }
//       }
//       if (found) break;
//     }
//     if (found) break;
//   }
//
//
//   PointType translate = pt0;
//
//   PointType x = pt1 - pt0;
//   x.normalize();
//
//   PointType y = pt2 - pt0;
//   y.normalize();
//
//   PointType z = viennagrid::cross_prod(x,y);
//   z.normalize();
//
//   y = viennagrid::cross_prod(x,z);
//   y.normalize();
//
//
//   Transformation from_2d;
//   Transformation to_2d;
//
//   from_2d.matrix.row_count = 3;
//   from_2d.matrix.column_count = 2;
//   from_2d.matrix.values.resize(6);
//   from_2d.matrix.values[0] = x[0];
//   from_2d.matrix.values[2] = x[1];
//   from_2d.matrix.values[4] = x[2];
//   from_2d.matrix.values[1] = y[0];
//   from_2d.matrix.values[3] = y[1];
//   from_2d.matrix.values[5] = y[2];
//   from_2d.translate = translate;
//
//   to_2d.matrix.row_count = 2;
//   to_2d.matrix.column_count = 3;
//   to_2d.matrix.values.resize(6);
//   to_2d.matrix.values[0] = x[0];
//   to_2d.matrix.values[1] = x[1];
//   to_2d.matrix.values[2] = x[2];
//   to_2d.matrix.values[3] = y[0];
//   to_2d.matrix.values[4] = y[1];
//   to_2d.matrix.values[5] = y[2];
//   to_2d.translate.resize(2);
//   to_2d.translate[0] = - viennagrid::inner_prod(x, translate);
//   to_2d.translate[1] = - viennagrid::inner_prod(y, translate);
//
//   return std::make_pair(to_2d, from_2d);
// }


std::pair<Transformation, Transformation> Transformation::make_facet_projection(viennagrid_plc plc, viennagrid_element_id facet_id, viennagrid_numeric eps)
{
  viennagrid_dimension geo_dim;
  viennagrid_plc_geometric_dimension_get(plc, &geo_dim);

  viennagrid_element_id * facet_vertices_begin;
  viennagrid_element_id * facet_vertices_end;
  viennagrid_plc_boundary_elements(plc, facet_id, 0, &facet_vertices_begin, &facet_vertices_end);
  viennagrid_int facet_vertex_count = facet_vertices_end-facet_vertices_begin;

  PointType centroid = PointType(geo_dim);
  PointType pt0;
  PointType pt1;

  viennagrid_element_id * vit0 = facet_vertices_begin;
  for (; vit0 != facet_vertices_end; ++vit0)
  {
    centroid += get_point(plc, *vit0);
  }
  centroid /= facet_vertex_count;


  vit0 = facet_vertices_begin;
  viennagrid_element_id * vit1 = vit0; ++vit1;

  PointType x;
  PointType y;

  bool found;
  for (; vit0 != facet_vertices_end; ++vit0)
  {
    viennagrid_element_id * vit1 = vit0; ++vit1;
    pt0 = get_point(plc, *vit0);
    for (; vit1 != facet_vertices_end; ++vit1)
    {
      pt1 = get_point(plc, *vit1);

      x = pt0 - centroid;
      y = pt1 - centroid;

      PointType c = viennagrid::cross_prod(x, y);
      if (viennagrid::norm_2(c) > eps)
      {
        found = true;
        break;
      }
    }
    if (found) break;
  }

  assert(found);


  PointType translate = centroid;

  x.normalize();
  y.normalize();

  PointType z = viennagrid::cross_prod(x,y);
  z.normalize();

  y = viennagrid::cross_prod(x,z);
  y.normalize();

//   std::cout << determinant(x,y,z) << std::endl;

  Transformation from_2d;
  Transformation to_2d;

  from_2d.matrix.row_count = 3;
  from_2d.matrix.column_count = 2;
  from_2d.matrix.values.resize(6);
  from_2d.matrix.values[0] = x[0];
  from_2d.matrix.values[2] = x[1];
  from_2d.matrix.values[4] = x[2];
  from_2d.matrix.values[1] = y[0];
  from_2d.matrix.values[3] = y[1];
  from_2d.matrix.values[5] = y[2];
  from_2d.translate = translate;

  to_2d.matrix.row_count = 2;
  to_2d.matrix.column_count = 3;
  to_2d.matrix.values.resize(6);
  to_2d.matrix.values[0] = x[0];
  to_2d.matrix.values[1] = x[1];
  to_2d.matrix.values[2] = x[2];
  to_2d.matrix.values[3] = y[0];
  to_2d.matrix.values[4] = y[1];
  to_2d.matrix.values[5] = y[2];
  to_2d.translate.resize(2);
  to_2d.translate[0] = - viennagrid::inner_prod(x, translate);
  to_2d.translate[1] = - viennagrid::inner_prod(y, translate);

  return std::make_pair(to_2d, from_2d);
}





MeshType plc_to_linemesh(viennagrid_plc plc)
{
  MeshType mesh;

  viennagrid_dimension geo_dim;
  viennagrid_plc_geometric_dimension_get(plc, &geo_dim);

  viennagrid_element_id vertices_begin;
  viennagrid_element_id vertices_end;
  viennagrid_plc_elements_get(plc, 0, &vertices_begin, &vertices_end);

  std::map<viennagrid_element_id, ElementType> vertex_map;
  for (viennagrid_element_id vid = vertices_begin; vid != vertices_end; ++vid)
  {
    viennagrid_numeric * pt;
    viennagrid_plc_vertex_coords_get(plc, vid, &pt);

    PointType mpt(geo_dim);
    std::copy(pt, pt+geo_dim, &mpt[0]);

    vertex_map[vid] = viennagrid::make_vertex(mesh, mpt);
  }

  viennagrid_element_id lines_begin;
  viennagrid_element_id lines_end;
  viennagrid_plc_elements_get(plc, 1, &lines_begin, &lines_end);
  for (viennagrid_element_id lid = lines_begin; lid != lines_end; ++lid)
  {
    viennagrid_element_id * line_vertices_begin;
    viennagrid_element_id * line_vertices_end;
    viennagrid_plc_boundary_elements(plc, lid, 0, &line_vertices_begin, &line_vertices_end);
    viennagrid::make_line(mesh, vertex_map[line_vertices_begin[0]], vertex_map[line_vertices_begin[1]]);
  }

  return mesh;
}






StatisticType boundary_min_angle_statistic(MeshType const & mesh)
{
  StatisticType statistic;
  statistic.set_histogram( StatisticType::histogram_type::make_uniform(0, M_PI, 18) );

  ElementRangeType vertices(mesh, 0);
  for (ElementRangeIterator vit = vertices.begin(); vit != vertices.end(); ++vit)
  {
    PointType pt = viennagrid::get_point(*vit);

    std::set<viennagrid_numeric> angles;

    NeighborRangeType neighbors(mesh, *vit, 1, 0);
    for (NeighborRangeIterator nvit = neighbors.begin(); nvit != neighbors.end(); ++nvit)
    {
      PointType npt = viennagrid::get_point(*nvit);

      PointType to_npt = npt-pt;
      to_npt.normalize();

      angles.insert( atan2(to_npt[1], to_npt[0]) );
    }

    std::set<viennagrid_numeric>::iterator prev_it = angles.begin();
    std::set<viennagrid_numeric>::iterator it = prev_it; ++it;

    for (; it != angles.end(); ++it, ++prev_it)
    {
      viennagrid_numeric angle = std::abs(*it - *prev_it);
      statistic(angle);
    }
  }

  return statistic;
}








QuantityField min_angles(MeshType const & mesh)
{
  QuantityField qf;
  int geo_dim = viennagrid::geometric_dimension(mesh);

  CellRangeType cells(mesh);
  for (CellRangeIterator cit = cells.begin(); cit != cells.end(); ++cit)
  {
    if (geo_dim == 2)
      qf.set(*cit, viennamesh::min_angle(*cit));
    else
      qf.set(*cit, viennamesh::min_dihedral_angle(*cit));
  }

  return qf;
}

QuantityField radius_edge_ratio(MeshType const & mesh)
{
  QuantityField qf;
  int geo_dim = viennagrid::geometric_dimension(mesh);

  CellRangeType cells(mesh);
  for (CellRangeIterator cit = cells.begin(); cit != cells.end(); ++cit)
  {
    qf.set(*cit, viennamesh::radius_edge_ratio(*cit));
  }

  return qf;
}

StatisticType statistics_from_quantity_field(MeshType const & mesh, QuantityField const & qf,
                                             viennagrid_numeric min, viennagrid_numeric max, std::size_t count)
{
  StatisticType statistic;
  statistic.set_histogram( StatisticType::histogram_type::make_uniform(min, max, count) );

  CellRangeType cells(mesh);
  for (CellRangeIterator cit = cells.begin(); cit != cells.end(); ++cit)
  {
    statistic( qf.get(*cit) );
  }

  return statistic;
}


int mesh_size(MeshType const & mesh)
{
  // Listing 2.1 (page 20)
  viennagrid_dimension geo_dim = viennagrid::geometric_dimension(mesh);
  int vertex_count = viennagrid::elements(mesh, 0).size();
  int cell_count = viennagrid::elements(mesh, 2).size(); // also is triangle count

  int total_size = 0;
  total_size += sizeof_int;                                                   // point_dimension
  total_size += sizeof_int;                                                   // vertex_count
  total_size += sizeof_ptr + sizeof_numeric * geo_dim * vertex_count;         // vertex_coords

  total_size += sizeof_int;                                                   // cell_count
  total_size += sizeof_ptr + sizeof_element_type * cell_count;                // cell_types
  total_size += sizeof_ptr + sizeof_index_type * (cell_count+1);              // cell_vertex_indices
  total_size += sizeof_ptr + sizeof_index_type * cell_count * 3;
  return total_size;
}

int mr_mesh_size(MeshType const & mesh)
{
  // Listing 2.2 (page 20)
  viennagrid_dimension geo_dim = viennagrid::geometric_dimension(mesh);
  int cell_count = viennagrid::elements(mesh, 2).size(); // also is triangle count

  int total_size = 0;
  total_size += sizeof_ptr + sizeof_region_id_type * cell_count;
  return total_size + mesh_size(mesh);                                        // cell_region
}



