#include "3d_geometry.hpp"



viennagrid_plc rotate_sketch_z(std::vector<PointType> const & pts,
                               std::vector< std::pair<int, int> > const & line_ids,
                               double angle, double eps)
{
  viennagrid_plc plc;
  viennagrid_plc_create(&plc);

  viennagrid_plc_geometric_dimension_set(plc, 3);

  std::size_t offset = pts.size();
  std::map<int, viennagrid_element_id> vertices;

  for (std::size_t i = 0; i != offset; ++i)
  {
    vertices[i] = make_vertex(plc, pts[i]);
  }

  for (std::size_t i = 0; i != offset; ++i)
  {
    PointType pt = pts[i];
    PointType rotated_pt = rotate_z(pt, angle);
    viennagrid_element_id rotated_vtx;
    if ( viennagrid::detail::is_equal(eps, pt, rotated_pt) )
      rotated_vtx = vertices[i];
    else
    {
      viennagrid_plc_vertex_create(plc, &rotated_pt[0], &rotated_vtx);
    }

    vertices[i+offset] = rotated_vtx;
  }

  std::vector<viennagrid_element_id> sketch_line_ids;
  std::vector<viennagrid_element_id> rotated_sketch_line_ids;

  for (std::size_t i = 0; i != line_ids.size(); ++i)
  {
    viennagrid_element_id line = make_line(plc, vertices[line_ids[i].first], vertices[line_ids[i].second]);
    viennagrid_element_id rotated_line;

    if ((vertices[line_ids[i].first] != vertices[line_ids[i].first+offset]) ||
        (vertices[line_ids[i].second] != vertices[line_ids[i].second+offset]))
    {
      rotated_line = make_line(plc, vertices[line_ids[i].first+offset], vertices[line_ids[i].second+offset]);
    }
    else
    {
      rotated_line = line;
    }

    sketch_line_ids.push_back(line);
    rotated_sketch_line_ids.push_back(rotated_line);
  }

  std::map<int, viennagrid_element_id> extruded_line_ids;

  for (std::size_t i = 0; i != offset; ++i)
  {
    if (vertices[i] != vertices[i+offset])
    {
//       viennagrid_element_id line;
//       viennagrid_plc_line_create(plc, vertices[i], vertices[i+offset], &line);
//       extruded_line_ids[i] = line;

      extruded_line_ids[i] = make_line(plc, vertices[i], vertices[i+offset]);
    }
  }

  viennagrid_plc_facet_create(plc, sketch_line_ids.size(), &sketch_line_ids[0], NULL);
  viennagrid_plc_facet_create(plc, rotated_sketch_line_ids.size(), &rotated_sketch_line_ids[0], NULL);


  std::vector< std::vector<int> > merged_line_ids;
  for (std::size_t i = 0; i != line_ids.size(); ++i)
  {
    bool found = false;

    for (std::size_t j = 0; j != merged_line_ids.size(); ++j)
    {
      PointType ld = pts[line_ids[merged_line_ids[j].front()].second] - pts[line_ids[merged_line_ids[j].front()].first];

      PointType d0 = pts[line_ids[i].first] - pts[line_ids[merged_line_ids[j].front()].first];
      PointType d1 = pts[line_ids[i].second] - pts[line_ids[merged_line_ids[j].front()].first];

      PointType c0 = viennagrid::cross_prod(ld, d0);
      PointType c1 = viennagrid::cross_prod(ld, d1);

      if ((viennagrid::norm_2(c0) < eps) && (viennagrid::norm_2(c1) < eps))
      {
//         std::cout << "Found two co-linear line_ids: " << pts[line_ids[merged_line_ids[j].front()].first] << " - " << pts[line_ids[merged_line_ids[j].front()].second] <<
//             "      with " <<  pts[line_ids[i].first] << " " << pts[line_ids[i].second] << std::endl;
        found = true;
        merged_line_ids[j].push_back(i);
        break;
      }
    }

    if (!found)
    {
      merged_line_ids.push_back( std::vector<int>(1,i) );
    }
  }

  for (std::size_t i = 0; i != merged_line_ids.size(); ++i)
  {
    std::vector<viennagrid_element_id> facet_line_ids;
    for (std::size_t j = 0; j != merged_line_ids[i].size(); ++j)
    {
      int line_id = merged_line_ids[i][j];

      facet_line_ids.push_back( sketch_line_ids[line_id] );
      if (sketch_line_ids[line_id] != rotated_sketch_line_ids[line_id])
        facet_line_ids.push_back(rotated_sketch_line_ids[line_id]);

      std::map<int, viennagrid_element_id>::iterator it = extruded_line_ids.find(line_ids[line_id].first);
      if (it != extruded_line_ids.end())
        facet_line_ids.push_back( (*it).second );

      it = extruded_line_ids.find(line_ids[line_id].second);
      if (it != extruded_line_ids.end())
        facet_line_ids.push_back( (*it).second );
    }

    if (facet_line_ids.size() > 2)
    {
      viennagrid_plc_facet_create(plc, facet_line_ids.size(), &facet_line_ids[0], NULL);
    }
  }


//   for (std::size_t i = 0; i != line_ids.size(); ++i)
//   {
//     if (sketch_line_ids[i] != rotated_sketch_line_ids[i])
//     {
//       std::vector<viennagrid_element_id> facet_line_ids;
//       facet_line_ids.push_back(sketch_line_ids[i]);
//       facet_line_ids.push_back(rotated_sketch_line_ids[i]);
//
//       std::map<int, viennagrid_element_id>::iterator it = extruded_line_ids.find(line_ids[i].first);
//       if (it != extruded_line_ids.end())
//         facet_line_ids.push_back( (*it).second );
//
//       it = extruded_line_ids.find(line_ids[i].second);
//       if (it != extruded_line_ids.end())
//         facet_line_ids.push_back( (*it).second );
//
//       viennagrid_plc_facet_create(plc, facet_line_ids.size(), &facet_line_ids[0], NULL);
//     }
//   }

  return plc;
}


Transformation find_transform(TriangleWrapper const & from,
                              TriangleWrapper const & to,
                              viennagrid_numeric eps,
                              bool swap_orientation)
{
  std::vector<PointType> from_pts;
  std::vector<PointType> to_pts;

  for (int i = 0; i != from.numberofpoints; ++i)
  {
    from_pts.push_back( viennagrid::make_point( from.pointlist[2*i+0], from.pointlist[2*i+1] ) );
  }

  for (int i = 0; i != from.numberofpoints; ++i)
  {
    to_pts.push_back( viennagrid::make_point( to.pointlist[2*i+0], to.pointlist[2*i+1] ) );
  }

  return Transformation::find_2d_rigid(from_pts, to_pts, eps, swap_orientation);
}



viennagrid_plc make_rotational_3d_geometry(TSVconfig const & tsv_cfg, viennagrid_numeric eps)
{
  std::vector<PointType> pts;
  pts.push_back( viennagrid::make_point(0,0,0) );
  pts.push_back( viennagrid::make_point(0,0,tsv_cfg.trunk_height()) );
  pts.push_back( viennagrid::make_point(0,0,tsv_cfg.trunk_height() + tsv_cfg.layer_thickness) );
  pts.push_back( viennagrid::make_point(tsv_cfg.radius,0,0) );
  pts.push_back( viennagrid::make_point(tsv_cfg.radius,0,tsv_cfg.trunk_height()) );
  pts.push_back( viennagrid::make_point(tsv_cfg.hole_radius,0,tsv_cfg.trunk_height() + tsv_cfg.layer_thickness) );
  pts.push_back( viennagrid::make_point(tsv_cfg.hole_radius + tsv_cfg.layer_thickness,0,tsv_cfg.trunk_height() + tsv_cfg.layer_thickness) );
  pts.push_back( viennagrid::make_point(tsv_cfg.radius,0,tsv_cfg.trunk_height() + tsv_cfg.layer_thickness) );
  pts.push_back( viennagrid::make_point(tsv_cfg.hole_radius + tsv_cfg.layer_thickness,0,tsv_cfg.total_height - tsv_cfg.layer_thickness) );
  pts.push_back( viennagrid::make_point(tsv_cfg.radius,0,tsv_cfg.total_height - tsv_cfg.layer_thickness) );
  pts.push_back( viennagrid::make_point(tsv_cfg.hole_radius,0,tsv_cfg.total_height) );
  pts.push_back( viennagrid::make_point(tsv_cfg.radius,0,tsv_cfg.total_height) );

  std::vector< std::pair<int, int> > line_ids;
  line_ids.push_back( std::make_pair(0,1) );
  line_ids.push_back( std::make_pair(1,2) );
  line_ids.push_back( std::make_pair(5,10) );
  line_ids.push_back( std::make_pair(6,8) );
  line_ids.push_back( std::make_pair(3,4) );
  line_ids.push_back( std::make_pair(4,7) );
  line_ids.push_back( std::make_pair(7,9) );
  line_ids.push_back( std::make_pair(9,11) );

  line_ids.push_back( std::make_pair(0,3) );
  line_ids.push_back( std::make_pair(1,4) );
  line_ids.push_back( std::make_pair(2,5) );
  line_ids.push_back( std::make_pair(6,7) );
  line_ids.push_back( std::make_pair(8,9) );
  line_ids.push_back( std::make_pair(10,11) );


  viennagrid_plc geometry = rotate_sketch_z(pts, line_ids, tsv_cfg.angle(), eps);

  viennamesh::seed_point_container seed_points;
  seed_points.push_back(viennamesh::seed_point(rotate_z(viennagrid::make_point(tsv_cfg.radius/2,0,tsv_cfg.trunk_height()/2), tsv_cfg.angle()/2), 0) );
  seed_points.push_back(viennamesh::seed_point(rotate_z(viennagrid::make_point((tsv_cfg.radius+tsv_cfg.hole_radius+tsv_cfg.layer_thickness)/2,0,tsv_cfg.trunk_height()+tsv_cfg.layer_thickness+tsv_cfg.body_height()/2), tsv_cfg.angle()/2), 0) );
  seed_points.push_back(viennamesh::seed_point(rotate_z(viennagrid::make_point(tsv_cfg.radius/2,0,tsv_cfg.trunk_height()+tsv_cfg.layer_thickness/2), tsv_cfg.angle()/2), 1) );
  add_seed_points(geometry, seed_points);

  return geometry;
}


viennagrid_plc make_rotational_3d_geometry_nonsym(TSVconfig const & tsv_cfg, viennagrid_numeric eps)
{
  std::vector<PointType> pts;
  pts.push_back( viennagrid::make_point(tsv_cfg.nonsym_inner_radius,0,0) );
  pts.push_back( viennagrid::make_point(tsv_cfg.nonsym_inner_radius,0,tsv_cfg.trunk_height()) );
  pts.push_back( viennagrid::make_point(tsv_cfg.nonsym_inner_radius,0,tsv_cfg.trunk_height() + tsv_cfg.layer_thickness) );
  pts.push_back( viennagrid::make_point(tsv_cfg.radius,0,0) );
  pts.push_back( viennagrid::make_point(tsv_cfg.radius,0,tsv_cfg.trunk_height()) );
  pts.push_back( viennagrid::make_point(tsv_cfg.hole_radius,0,tsv_cfg.trunk_height() + tsv_cfg.layer_thickness) );
  pts.push_back( viennagrid::make_point(tsv_cfg.hole_radius + tsv_cfg.layer_thickness,0,tsv_cfg.trunk_height() + tsv_cfg.layer_thickness) );
  pts.push_back( viennagrid::make_point(tsv_cfg.radius,0,tsv_cfg.trunk_height() + tsv_cfg.layer_thickness) );
  pts.push_back( viennagrid::make_point(tsv_cfg.hole_radius + tsv_cfg.layer_thickness,0,tsv_cfg.total_height - tsv_cfg.layer_thickness) );
  pts.push_back( viennagrid::make_point(tsv_cfg.radius,0,tsv_cfg.total_height - tsv_cfg.layer_thickness) );
  pts.push_back( viennagrid::make_point(tsv_cfg.hole_radius,0,tsv_cfg.total_height) );
  pts.push_back( viennagrid::make_point(tsv_cfg.radius,0,tsv_cfg.total_height) );

  std::vector< std::pair<int, int> > line_ids;
  line_ids.push_back( std::make_pair(0,1) );
  line_ids.push_back( std::make_pair(1,2) );
  line_ids.push_back( std::make_pair(5,10) );
  line_ids.push_back( std::make_pair(6,8) );
  line_ids.push_back( std::make_pair(3,4) );
  line_ids.push_back( std::make_pair(4,7) );
  line_ids.push_back( std::make_pair(7,9) );
  line_ids.push_back( std::make_pair(9,11) );

  line_ids.push_back( std::make_pair(0,3) );
  line_ids.push_back( std::make_pair(1,4) );
  line_ids.push_back( std::make_pair(2,5) );
  line_ids.push_back( std::make_pair(6,7) );
  line_ids.push_back( std::make_pair(8,9) );
  line_ids.push_back( std::make_pair(10,11) );


  viennagrid_plc geometry = rotate_sketch_z(pts, line_ids, tsv_cfg.angle(), eps);

  viennamesh::seed_point_container seed_points;
  seed_points.push_back(viennamesh::seed_point(rotate_z(viennagrid::make_point((tsv_cfg.nonsym_inner_radius+tsv_cfg.radius)/2,0,tsv_cfg.trunk_height()/2), tsv_cfg.angle()/2), 0) );
  seed_points.push_back(viennamesh::seed_point(rotate_z(viennagrid::make_point((tsv_cfg.radius+tsv_cfg.hole_radius+tsv_cfg.layer_thickness)/2,0,tsv_cfg.trunk_height()+tsv_cfg.layer_thickness+tsv_cfg.body_height()/2), tsv_cfg.angle()/2), 0) );
  seed_points.push_back(viennamesh::seed_point(rotate_z(viennagrid::make_point(tsv_cfg.radius/2,0,tsv_cfg.trunk_height()+tsv_cfg.layer_thickness/2), tsv_cfg.angle()/2), 1) );
  add_seed_points(geometry, seed_points);

  return geometry;
}

viennagrid_plc make_rotational_3d_geometry_nonsym_middle(TSVconfig const & tsv_cfg)
{
  viennagrid_plc nonsym_plc = make_plc(3);

  std::vector<viennagrid_element_id> vertices(3*tsv_cfg.rotational_symmetry_order);

  for (int i = 0; i != tsv_cfg.rotational_symmetry_order; ++i)
  {
    vertices[3*i+0] = make_vertex(nonsym_plc, rotate_z(viennagrid::make_point(tsv_cfg.nonsym_inner_radius,0,0), tsv_cfg.angle()*i));
    vertices[3*i+1] = make_vertex(nonsym_plc, rotate_z(viennagrid::make_point(tsv_cfg.nonsym_inner_radius,0,tsv_cfg.trunk_height()), tsv_cfg.angle()*i));
    vertices[3*i+2] = make_vertex(nonsym_plc, rotate_z(viennagrid::make_point(tsv_cfg.nonsym_inner_radius,0,tsv_cfg.trunk_height() + tsv_cfg.layer_thickness), tsv_cfg.angle()*i));
  }
  std::vector<viennagrid_element_id> line_ids(5*tsv_cfg.rotational_symmetry_order);
  for (int i = 0; i != tsv_cfg.rotational_symmetry_order; ++i)
  {
    line_ids[5*i+0] = make_line(nonsym_plc, vertices[3*i+0], vertices[3*i+1]);
    line_ids[5*i+1] = make_line(nonsym_plc, vertices[3*i+1], vertices[3*i+2]);

    line_ids[5*i+2] = make_line(nonsym_plc, vertices[3*i+0], vertices[3*((i+1)%tsv_cfg.rotational_symmetry_order)+0]);
    line_ids[5*i+3] = make_line(nonsym_plc, vertices[3*i+1], vertices[3*((i+1)%tsv_cfg.rotational_symmetry_order)+1]);
    line_ids[5*i+4] = make_line(nonsym_plc, vertices[3*i+2], vertices[3*((i+1)%tsv_cfg.rotational_symmetry_order)+2]);
  }

  for (int i = 0; i != tsv_cfg.rotational_symmetry_order; ++i)
  {
    std::vector<viennagrid_element_id> facet_line_ids;
    facet_line_ids.push_back(line_ids[5*i+0]);
    facet_line_ids.push_back(line_ids[5*i+1]);
    facet_line_ids.push_back(line_ids[5*i+2]);
    facet_line_ids.push_back(line_ids[5*i+3]);
    facet_line_ids.push_back(line_ids[5*i+4]);
    facet_line_ids.push_back(line_ids[5*((i+1)%tsv_cfg.rotational_symmetry_order)+0]);
    facet_line_ids.push_back(line_ids[5*((i+1)%tsv_cfg.rotational_symmetry_order)+1]);
    viennagrid_plc_facet_create(nonsym_plc, facet_line_ids.size(), &facet_line_ids[0], NULL);
  }

  {
    std::vector<viennagrid_element_id> facet_line_ids;
    for (int i = 0; i != tsv_cfg.rotational_symmetry_order; ++i)
      facet_line_ids.push_back(line_ids[5*i+2]);
    viennagrid_plc_facet_create(nonsym_plc, facet_line_ids.size(), &facet_line_ids[0], NULL);
  }

  {
    std::vector<viennagrid_element_id> facet_line_ids;
    for (int i = 0; i != tsv_cfg.rotational_symmetry_order; ++i)
      facet_line_ids.push_back(line_ids[5*i+3]);
    viennagrid_plc_facet_create(nonsym_plc, facet_line_ids.size(), &facet_line_ids[0], NULL);
  }

  {
    std::vector<viennagrid_element_id> facet_line_ids;
    for (int i = 0; i != tsv_cfg.rotational_symmetry_order; ++i)
      facet_line_ids.push_back(line_ids[5*i+4]);
    viennagrid_plc_facet_create(nonsym_plc, facet_line_ids.size(), &facet_line_ids[0], NULL);
  }

  viennamesh::seed_point_container seed_points;
  seed_points.push_back(viennamesh::seed_point(rotate_z(viennagrid::make_point(0,0,tsv_cfg.trunk_height()/2), tsv_cfg.angle()/2), 0) );
  seed_points.push_back(viennamesh::seed_point(rotate_z(viennagrid::make_point(0,0,tsv_cfg.trunk_height()+tsv_cfg.layer_thickness/2), tsv_cfg.angle()/2), 1) );
  add_seed_points(nonsym_plc, seed_points);

  return nonsym_plc;
}



viennagrid_plc make_multi_tsv_geometry(TSVconfig const & tsv_cfg, std::vector<PointType> const & centers, viennagrid_numeric eps)
{
  int tsv_count = centers.size();
  viennagrid_numeric distance = tsv_cfg.distance;
  viennagrid_numeric trunk_size = tsv_cfg.trunk_size;


// //   PointType centers[tsv_count];
//
//   centers[0] = viennagrid::make_point(0,0,0);
//   centers[1] = viennagrid::make_point(0,1,0) * distance;
//   centers[2] = viennagrid::make_point(std::cos(M_PI/6),std::sin(M_PI/6),0) * distance;
//   centers[3] = viennagrid::make_point(std::cos(M_PI/6),-std::sin(M_PI/6),0) * distance;
//   centers[4] = viennagrid::make_point(0,-1,0) * distance;
//   centers[5] = viennagrid::make_point(-std::cos(M_PI/6),-std::sin(M_PI/6),0) * distance;
//   centers[6] = viennagrid::make_point(-std::cos(M_PI/6),std::sin(M_PI/6),0) * distance;


  std::vector<viennagrid_numeric> heights(5);
  heights[0] = 0;
  heights[1] = tsv_cfg.trunk_height();
  heights[2] = tsv_cfg.trunk_height()+tsv_cfg.layer_thickness;
  heights[3] = tsv_cfg.total_height-tsv_cfg.layer_thickness;
  heights[4] = tsv_cfg.total_height;

  int hc = heights.size();


  viennagrid_plc plc = make_plc(3);

  std::vector< std::vector<viennagrid_element_id> > tsv_vertices(tsv_count);
  std::vector< std::vector<viennagrid_element_id> > tsv_line_ids(tsv_count);
  for (int tsv = 0; tsv != tsv_count; ++tsv)
  {
    std::vector<viennagrid_element_id> & vertices = tsv_vertices[tsv];
    std::vector<viennagrid_element_id> & line_ids = tsv_line_ids[tsv];
    PointType const & center = centers[tsv];

    vertices.resize(hc*tsv_cfg.rotational_symmetry_order);
    for (int i = 0; i != tsv_cfg.rotational_symmetry_order; ++i)
    {
      for (int j = 0; j != hc; ++j)
      {
        vertices[5*i+j] = make_vertex(plc, center + rotate_z(viennagrid::make_point(tsv_cfg.radius,0,heights[j]), tsv_cfg.angle()*i));
      }
    }

    line_ids.resize((2*hc-1)*tsv_cfg.rotational_symmetry_order);
    for (int i = 0; i != tsv_cfg.rotational_symmetry_order; ++i)
    {
      for (int j = 0; j != hc-1; ++j)
        line_ids[9*i+j] = make_line(plc, vertices[5*i+j], vertices[5*i+j+1]);

      for (int j = 0; j != hc; ++j)
        line_ids[9*i+hc-1+j] = make_line(plc, vertices[5*i+j], vertices[5*((i+1)%tsv_cfg.rotational_symmetry_order)+j]);
    }

    for (int i = 0; i != tsv_cfg.rotational_symmetry_order; ++i)
    {
      std::vector<viennagrid_element_id> facet_line_ids;
      for (int j = 0; j != 2*hc-1; ++j)
        facet_line_ids.push_back(line_ids[9*i+j]);

      for (int j = 0; j != hc-1; ++j)
        facet_line_ids.push_back(line_ids[9*((i+1)%tsv_cfg.rotational_symmetry_order)+j]);

      viennagrid_plc_facet_create(plc, facet_line_ids.size(), &facet_line_ids[0], NULL);
    }
  }






  std::vector<viennagrid_element_id> border_vertices(hc*4);
  for (int j = 0; j != hc; ++j)
  {
    border_vertices[     j] = make_vertex(plc, viennagrid::make_point(-trunk_size,-trunk_size,heights[j]));
    border_vertices[  hc+j] = make_vertex(plc, viennagrid::make_point( trunk_size,-trunk_size,heights[j]));
    border_vertices[2*hc+j] = make_vertex(plc, viennagrid::make_point( trunk_size, trunk_size,heights[j]));
    border_vertices[3*hc+j] = make_vertex(plc, viennagrid::make_point(-trunk_size, trunk_size,heights[j]));
  }

  std::vector<viennagrid_element_id> border_line_ids(4*hc + 4*(hc-1));

  for (int j = 0; j != hc-1; ++j)
  {
    for (int i = 0; i != 4; ++i)
    {
      border_line_ids[i*(2*hc-1) + j] = make_line(plc, border_vertices[i*hc+j], border_vertices[i*hc+j+1]);
      border_line_ids[i*(2*hc-1) + hc-1 + j] = make_line(plc, border_vertices[hc*i+j], border_vertices[hc*((i+1)%4)+j]);
    }

//       border_line_ids[             j] = make_line(plc, border_vertices[     j], border_vertices[     j+1]);
//     border_line_ids[   2*hc-1  + j] = make_line(plc, border_vertices[  hc+j], border_vertices[  hc+j+1]);
//     border_line_ids[2*(2*hc-1) + j] = make_line(plc, border_vertices[2*hc+j], border_vertices[2*hc+j+1]);
//     border_line_ids[3*(2*hc-1) + j] = make_line(plc, border_vertices[3*hc+j], border_vertices[3*hc+j+1]);
  }

  for (int j = 0; j != hc; ++j)
  {
    for (int i = 0; i != 4; ++i)
      border_line_ids[i*(2*hc-1) + hc-1 + j] = make_line(plc, border_vertices[hc*i+j], border_vertices[hc*((i+1)%4)+j]);
  }

  for (int j = 0; j != hc; ++j)
  {
    std::vector<viennagrid_element_id> facet_line_ids;

    for (int tsv = 0; tsv != tsv_count; ++tsv)
    {
      for (int i = 0; i != tsv_cfg.rotational_symmetry_order; ++i)
        facet_line_ids.push_back( tsv_line_ids[tsv][9*i+hc-1+j] );
    }

    facet_line_ids.push_back(border_line_ids[0*(2*hc-1) + hc-1 + j]);
    facet_line_ids.push_back(border_line_ids[1*(2*hc-1) + hc-1 + j]);
    facet_line_ids.push_back(border_line_ids[2*(2*hc-1) + hc-1 + j]);
    facet_line_ids.push_back(border_line_ids[3*(2*hc-1) + hc-1 + j]);

    viennagrid_element_id facet;
    viennagrid_plc_facet_create(plc, facet_line_ids.size(), &facet_line_ids[0], &facet);

    for (int tsv = 0; tsv != tsv_count; ++tsv)
    {
      PointType hole_point = centers[tsv];
      hole_point[2] = heights[j];
      viennagrid_plc_facet_hole_point_add(plc, facet, &hole_point[0]);
    }
  }

  for (int i = 0; i != 4; ++i)
  {
    std::vector<viennagrid_element_id> facet_line_ids;

    for (int j = 0; j != hc-1; ++j)
    {
      facet_line_ids.push_back( border_line_ids[        i*(2*hc-1) + j] );
      facet_line_ids.push_back( border_line_ids[((i+1)%4)*(2*hc-1) + j] );
    }

    for (int j = 0; j != hc; ++j)
    {
      facet_line_ids.push_back( border_line_ids[i*(2*hc-1) + hc-1 + j] );
    }

    viennagrid_plc_facet_create(plc, facet_line_ids.size(), &facet_line_ids[0], NULL);
  }


  viennamesh::seed_point_container seed_points;
  seed_points.push_back(viennamesh::seed_point(viennagrid::make_point(-trunk_size+1,-trunk_size+1,(heights[0]+heights[1])/2), 0));
  seed_points.push_back(viennamesh::seed_point(viennagrid::make_point(-trunk_size+1,-trunk_size+1,(heights[1]+heights[2])/2), 1));
  seed_points.push_back(viennamesh::seed_point(viennagrid::make_point(-trunk_size+1,-trunk_size+1,(heights[2]+heights[3])/2), 0));
  seed_points.push_back(viennamesh::seed_point(viennagrid::make_point(-trunk_size+1,-trunk_size+1,(heights[3]+heights[4])/2), 1));
  add_seed_points(plc, seed_points);


//   write_mesh(plc_to_linemesh(plc), "multi_tsv");
//
//   viennamesh::algorithm_handle surface_mesher = context.make_algorithm("triangle_make_hull");
//   surface_mesher.set_input( "geometry", plc );
//   surface_mesher.run();

//   viennamesh::algorithm_handle mesher = context.make_algorithm("tetgen_make_mesh");
//   mesher.link_input("geometry", surface_mesher, "mesh");
//   mesher.run();
//
//   write_mesh(surface_mesher.get_output<viennagrid_mesh>("mesh")(), "multi_tsv_surface");

  return plc;
}



// viennagrid_plc make_beam(BridgeConfig const & cfg, viennagrid_numeric eps)
// {
//   viennagrid_plc plc = make_plc(3);
//
//   std::vector<PointType> pts;
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2, -cfg.beam_size_y/2, -cfg.beam_length/2-cfg.beam_trunk_size) );
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2,  cfg.beam_size_y/2, -cfg.beam_length/2-cfg.beam_trunk_size) );
//   pts.push_back( viennagrid::make_point( cfg.beam_size_x/2, -cfg.beam_size_y/2, -cfg.beam_length/2-cfg.beam_trunk_size) );
//   pts.push_back( viennagrid::make_point( cfg.beam_size_x/2,  cfg.beam_size_y/2, -cfg.beam_length/2-cfg.beam_trunk_size) );
//
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2, -cfg.beam_size_y/2, -cfg.beam_length/2) );
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2,  cfg.beam_size_y/2, -cfg.beam_length/2) );
//   pts.push_back( viennagrid::make_point( cfg.beam_size_x/2, -cfg.beam_size_y/2, -cfg.beam_length/2) );
//   pts.push_back( viennagrid::make_point( cfg.beam_size_x/2,  cfg.beam_size_y/2, -cfg.beam_length/2) );
//
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2, -cfg.beam_size_y/2+cfg.beam_offset, -cfg.beam_length/2) );
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2,  cfg.beam_size_y/2-cfg.beam_offset, -cfg.beam_length/2) );
//   pts.push_back( viennagrid::make_point( cfg.beam_size_x/2, -cfg.beam_size_y/2+cfg.beam_offset, -cfg.beam_length/2) );
//   pts.push_back( viennagrid::make_point( cfg.beam_size_x/2,  cfg.beam_size_y/2-cfg.beam_offset, -cfg.beam_length/2) );
//
//   pts.push_back( viennagrid::make_point(-cfg.beam_inner_thickness/2, -cfg.beam_inner_height/2, -cfg.beam_length/2) );
//   pts.push_back( viennagrid::make_point(-cfg.beam_inner_thickness/2,  cfg.beam_inner_height/2, -cfg.beam_length/2) );
//   pts.push_back( viennagrid::make_point( cfg.beam_inner_thickness/2, -cfg.beam_inner_height/2, -cfg.beam_length/2) );
//   pts.push_back( viennagrid::make_point( cfg.beam_inner_thickness/2,  cfg.beam_inner_height/2, -cfg.beam_length/2) );
//
//
//
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2, -cfg.beam_size_y/2, cfg.beam_length/2+cfg.beam_trunk_size) );
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2,  cfg.beam_size_y/2, cfg.beam_length/2+cfg.beam_trunk_size) );
//   pts.push_back( viennagrid::make_point( cfg.beam_size_x/2, -cfg.beam_size_y/2, cfg.beam_length/2+cfg.beam_trunk_size) );
//   pts.push_back( viennagrid::make_point( cfg.beam_size_x/2,  cfg.beam_size_y/2, cfg.beam_length/2+cfg.beam_trunk_size) );
//
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2, -cfg.beam_size_y/2, cfg.beam_length/2) );
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2,  cfg.beam_size_y/2, cfg.beam_length/2) );
//   pts.push_back( viennagrid::make_point( cfg.beam_size_x/2, -cfg.beam_size_y/2, cfg.beam_length/2) );
//   pts.push_back( viennagrid::make_point( cfg.beam_size_x/2,  cfg.beam_size_y/2, cfg.beam_length/2) );
//
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2, -cfg.beam_size_y/2+cfg.beam_offset, cfg.beam_length/2) );
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2,  cfg.beam_size_y/2-cfg.beam_offset, cfg.beam_length/2) );
//   pts.push_back( viennagrid::make_point( cfg.beam_size_x/2, -cfg.beam_size_y/2+cfg.beam_offset, cfg.beam_length/2) );
//   pts.push_back( viennagrid::make_point( cfg.beam_size_x/2,  cfg.beam_size_y/2-cfg.beam_offset, cfg.beam_length/2) );
//
//   pts.push_back( viennagrid::make_point(-cfg.beam_inner_thickness/2, -cfg.beam_inner_height/2, cfg.beam_length/2) );
//   pts.push_back( viennagrid::make_point(-cfg.beam_inner_thickness/2,  cfg.beam_inner_height/2, cfg.beam_length/2) );
//   pts.push_back( viennagrid::make_point( cfg.beam_inner_thickness/2, -cfg.beam_inner_height/2, cfg.beam_length/2) );
//   pts.push_back( viennagrid::make_point( cfg.beam_inner_thickness/2,  cfg.beam_inner_height/2, cfg.beam_length/2) );
//
//   for (std::size_t i = 0; i != pts.size(); ++i)
//     make_vertex(plc, pts[i]);
//
//   std::vector< std::pair<int, int> > line_ids;
//   line_ids.push_back( std::make_pair(0,1) );
//   line_ids.push_back( std::make_pair(0,2) );
//   line_ids.push_back( std::make_pair(1,3) );
//   line_ids.push_back( std::make_pair(2,3) );
//
//   line_ids.push_back( std::make_pair(0,4) );
//   line_ids.push_back( std::make_pair(1,5) );
//   line_ids.push_back( std::make_pair(2,6) );
//   line_ids.push_back( std::make_pair(3,7) );
//
//   line_ids.push_back( std::make_pair(4,6) );
//   line_ids.push_back( std::make_pair(5,7) );
//
//   line_ids.push_back( std::make_pair(4,8) );
//   line_ids.push_back( std::make_pair(8,9) );
//   line_ids.push_back( std::make_pair(9,5) );
//   line_ids.push_back( std::make_pair(8,12) );
//   line_ids.push_back( std::make_pair(12,13) );
//   line_ids.push_back( std::make_pair(13,9) );
//
//   line_ids.push_back( std::make_pair(6,10) );
//   line_ids.push_back( std::make_pair(10,11) );
//   line_ids.push_back( std::make_pair(11,7) );
//   line_ids.push_back( std::make_pair(10,14) );
//   line_ids.push_back( std::make_pair(14,15) );
//   line_ids.push_back( std::make_pair(15,11) );
//
//
//   line_ids.push_back( std::make_pair(16+0,16+1) );
//   line_ids.push_back( std::make_pair(16+0,16+2) );
//   line_ids.push_back( std::make_pair(16+1,16+3) );
//   line_ids.push_back( std::make_pair(16+2,16+3) );
//
//   line_ids.push_back( std::make_pair(16+0,16+4) );
//   line_ids.push_back( std::make_pair(16+1,16+5) );
//   line_ids.push_back( std::make_pair(16+2,16+6) );
//   line_ids.push_back( std::make_pair(16+3,16+7) );
//
//   line_ids.push_back( std::make_pair(16+4,16+6) );
//   line_ids.push_back( std::make_pair(16+5,16+7) );
//
//   line_ids.push_back( std::make_pair(16+4,16+8) );
//   line_ids.push_back( std::make_pair(16+8,16+9) );
//   line_ids.push_back( std::make_pair(16+9,16+5) );
//   line_ids.push_back( std::make_pair(16+8,16+12) );
//   line_ids.push_back( std::make_pair(16+12,16+13) );
//   line_ids.push_back( std::make_pair(16+13,16+9) );
//
//   line_ids.push_back( std::make_pair(16+6,16+10) );
//   line_ids.push_back( std::make_pair(16+10,16+11) );
//   line_ids.push_back( std::make_pair(16+11,16+7) );
//   line_ids.push_back( std::make_pair(16+10,16+14) );
//   line_ids.push_back( std::make_pair(16+14,16+15) );
//   line_ids.push_back( std::make_pair(16+15,16+11) );
//
//
//   line_ids.push_back( std::make_pair(4,16+4) );
//   line_ids.push_back( std::make_pair(5,16+5) );
//   line_ids.push_back( std::make_pair(6,16+6) );
//   line_ids.push_back( std::make_pair(7,16+7) );
//   line_ids.push_back( std::make_pair(8,16+8) );
//   line_ids.push_back( std::make_pair(9,16+9) );
//   line_ids.push_back( std::make_pair(10,16+10) );
//   line_ids.push_back( std::make_pair(11,16+11) );
//   line_ids.push_back( std::make_pair(12,16+12) );
//   line_ids.push_back( std::make_pair(13,16+13) );
//   line_ids.push_back( std::make_pair(14,16+14) );
//   line_ids.push_back( std::make_pair(15,16+15) );
//
//   std::vector<viennagrid_element_id> lines(line_ids.size());
//   for (std::size_t i = 0; i != line_ids.size(); ++i)
//     lines[i] = make_line(plc, line_ids[i].first, line_ids[i].second);
//
//
//   for (int i = 0; i != 2; ++i)
//   {
//     int lo = 22*i;
//     {
//       viennagrid_element_id facet_line_ids[] = {lo+0,lo+1,lo+2,lo+3};
//       for (int i = 0; i != 4; ++i) facet_line_ids[i] = lines[facet_line_ids[i]];
//       viennagrid_plc_facet_create(plc, 4, facet_line_ids, NULL);
//     }
//
//     {
//       viennagrid_element_id facet_line_ids[] = {lo+0,lo+4,lo+10,lo+11,lo+12,lo+5};
//       for (int i = 0; i != 6; ++i) facet_line_ids[i] = lines[facet_line_ids[i]];
//       viennagrid_plc_facet_create(plc, 6, facet_line_ids, NULL);
//     }
//     {
//       viennagrid_element_id facet_line_ids[] = {lo+5,lo+9,lo+2,lo+7};
//       for (int i = 0; i != 4; ++i) facet_line_ids[i] = lines[facet_line_ids[i]];
//       viennagrid_plc_facet_create(plc, 4, facet_line_ids, NULL);
//     }
//     {
//       viennagrid_element_id facet_line_ids[] = {lo+3,lo+6,lo+7,lo+16,lo+17,lo+18};
//       for (int i = 0; i != 6; ++i) facet_line_ids[i] = lines[facet_line_ids[i]];
//       viennagrid_plc_facet_create(plc, 6, facet_line_ids, NULL);
//     }
//     {
//       viennagrid_element_id facet_line_ids[] = {lo+4,lo+8,lo+6,lo+1};
//       for (int i = 0; i != 4; ++i) facet_line_ids[i] = lines[facet_line_ids[i]];
//       viennagrid_plc_facet_create(plc, 4, facet_line_ids, NULL);
//     }
//
//     {
//       viennagrid_element_id facet_line_ids[] = {lo+8,lo+9,lo+10,lo+11,lo+12,lo+13,lo+14,lo+15,lo+16,lo+17,lo+18,lo+19,lo+20,lo+21};
//       for (int i = 0; i != 14; ++i) facet_line_ids[i] = lines[facet_line_ids[i]];
//       viennagrid_plc_facet_create(plc, 14, facet_line_ids, NULL);
//     }
//   }
//
//   {
//     viennagrid_element_id facet_line_ids[] = {8,22+8,40+4,40+6};
//     for (int i = 0; i != 4; ++i) facet_line_ids[i] = lines[facet_line_ids[i]];
//     viennagrid_plc_facet_create(plc, 4, facet_line_ids, NULL);
//   }
//   {
//     viennagrid_element_id facet_line_ids[] = {9,22+9,40+5,40+7};
//     for (int i = 0; i != 4; ++i) facet_line_ids[i] = lines[facet_line_ids[i]];
//     viennagrid_plc_facet_create(plc, 4, facet_line_ids, NULL);
//   }
//
//   {
//     viennagrid_element_id facet_line_ids[] = {10,22+10,40+4,40+8};
//     for (int i = 0; i != 4; ++i) facet_line_ids[i] = lines[facet_line_ids[i]];
//     viennagrid_plc_facet_create(plc, 4, facet_line_ids, NULL);
//   }
//   {
//     viennagrid_element_id facet_line_ids[] = {12,22+12,40+5,40+9};
//     for (int i = 0; i != 4; ++i) facet_line_ids[i] = lines[facet_line_ids[i]];
//     viennagrid_plc_facet_create(plc, 4, facet_line_ids, NULL);
//   }
//   {
//     viennagrid_element_id facet_line_ids[] = {13,22+13,40+8,40+12};
//     for (int i = 0; i != 4; ++i) facet_line_ids[i] = lines[facet_line_ids[i]];
//     viennagrid_plc_facet_create(plc, 4, facet_line_ids, NULL);
//   }
//   {
//     viennagrid_element_id facet_line_ids[] = {14,22+14,40+12,40+13};
//     for (int i = 0; i != 4; ++i) facet_line_ids[i] = lines[facet_line_ids[i]];
//     viennagrid_plc_facet_create(plc, 4, facet_line_ids, NULL);
//   }
//   {
//     viennagrid_element_id facet_line_ids[] = {15,22+15,40+9,40+13};
//     for (int i = 0; i != 4; ++i) facet_line_ids[i] = lines[facet_line_ids[i]];
//     viennagrid_plc_facet_create(plc, 4, facet_line_ids, NULL);
//   }
//
//   {
//     viennagrid_element_id facet_line_ids[] = {16,22+16,40+6,40+10};
//     for (int i = 0; i != 4; ++i) facet_line_ids[i] = lines[facet_line_ids[i]];
//     viennagrid_plc_facet_create(plc, 4, facet_line_ids, NULL);
//   }
//   {
//     viennagrid_element_id facet_line_ids[] = {18,22+18,40+7,40+11};
//     for (int i = 0; i != 4; ++i) facet_line_ids[i] = lines[facet_line_ids[i]];
//     viennagrid_plc_facet_create(plc, 4, facet_line_ids, NULL);
//   }
//   {
//     viennagrid_element_id facet_line_ids[] = {19,22+19,40+10,40+14};
//     for (int i = 0; i != 4; ++i) facet_line_ids[i] = lines[facet_line_ids[i]];
//     viennagrid_plc_facet_create(plc, 4, facet_line_ids, NULL);
//   }
//   {
//     viennagrid_element_id facet_line_ids[] = {20,22+20,40+14,40+15};
//     for (int i = 0; i != 4; ++i) facet_line_ids[i] = lines[facet_line_ids[i]];
//     viennagrid_plc_facet_create(plc, 4, facet_line_ids, NULL);
//   }
//   {
//     viennagrid_element_id facet_line_ids[] = {21,22+21,40+15,40+11};
//     for (int i = 0; i != 4; ++i) facet_line_ids[i] = lines[facet_line_ids[i]];
//     viennagrid_plc_facet_create(plc, 4, facet_line_ids, NULL);
//   }
//
//   return plc;
// }



// viennagrid_plc make_beam_connector(BridgeConfig const & cfg, viennagrid_numeric eps)
// {
//   viennagrid_plc plc = make_plc(3);
//
//   std::vector<PointType> pts;
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2, 0, -cfg.beam_size_y*sqrt(3)/2) );
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2, 0,  cfg.beam_size_y*sqrt(3)/2) );
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2, cfg.beam_size_y,  cfg.beam_size_y*sqrt(3)/2) );
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2, cfg.beam_size_y*3/2, 0) );
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2, cfg.beam_size_y, -cfg.beam_size_y*sqrt(3)/2) );
//
//   pts.push_back( viennagrid::make_point( cfg.beam_size_x/2, 0, -cfg.beam_size_y*sqrt(3)/2) );
//   pts.push_back( viennagrid::make_point( cfg.beam_size_x/2, 0,  cfg.beam_size_y*sqrt(3)/2) );
//   pts.push_back( viennagrid::make_point( cfg.beam_size_x/2, cfg.beam_size_y,  cfg.beam_size_y*sqrt(3)/2) );
//   pts.push_back( viennagrid::make_point( cfg.beam_size_x/2, cfg.beam_size_y*3/2, 0) );
//   pts.push_back( viennagrid::make_point( cfg.beam_size_x/2, cfg.beam_size_y, -cfg.beam_size_y*sqrt(3)/2) );
//
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2, 0, -cfg.beam_size_x/2) );
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2, 0, +cfg.beam_size_x/2) );
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2, cfg.beam_size_y, -cfg.beam_size_x/2) );
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2, cfg.beam_size_y, +cfg.beam_size_x/2) );
//
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2-cfg.beam_trunk_size, 0, -cfg.beam_size_x/2) );
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2-cfg.beam_trunk_size, 0, +cfg.beam_size_x/2) );
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2-cfg.beam_trunk_size, cfg.beam_size_y, -cfg.beam_size_x/2) );
//   pts.push_back( viennagrid::make_point(-cfg.beam_size_x/2-cfg.beam_trunk_size, cfg.beam_size_y, +cfg.beam_size_x/2) );
//
//   for (std::size_t i = 0; i != pts.size(); ++i)
//     make_vertex(plc, pts[i]);
//
//   std::vector<viennagrid_element_id> lines;
//   lines.push_back( make_line(plc, 0, 10) );   // 0
//   lines.push_back( make_line(plc, 10, 10) );
//   lines.push_back( make_line(plc, 11, 1) );
//   lines.push_back( make_line(plc, 1, 2) );
//   lines.push_back( make_line(plc, 2, 3) );
//   lines.push_back( make_line(plc, 3, 4) );
//   lines.push_back( make_line(plc, 4, 0) );
//
//   lines.push_back( make_line(plc, 5, 6) );    // 7
//   lines.push_back( make_line(plc, 6, 7) );
//   lines.push_back( make_line(plc, 7, 8) );
//   lines.push_back( make_line(plc, 8, 9) );
//   lines.push_back( make_line(plc, 9, 5) );
//
//   lines.push_back( make_line(plc, 0, 5) );    // 12
//   lines.push_back( make_line(plc, 1, 6) );
//   lines.push_back( make_line(plc, 2, 7) );
//   lines.push_back( make_line(plc, 3, 8) );
//   lines.push_back( make_line(plc, 4, 9) );
//
//   lines.push_back( make_line(plc, 10, 12) );  // 17
//   lines.push_back( make_line(plc, 12, 13) );
//   lines.push_back( make_line(plc, 11, 13) );
//
//   lines.push_back( make_line(plc, 14, 15) );  // 20
//   lines.push_back( make_line(plc, 14, 16) );
//   lines.push_back( make_line(plc, 16, 17) );
//   lines.push_back( make_line(plc, 15, 17) );
//
//   lines.push_back( make_line(plc, 10, 14) );  // 24
//   lines.push_back( make_line(plc, 11, 15) );
//   lines.push_back( make_line(plc, 12, 16) );
//   lines.push_back( make_line(plc, 13, 17) );
//
//   {
//     viennagrid_element_id facet_line_ids[] = {lines[20], lines[21], lines[22], lines[23]};
//     viennagrid_plc_facet_create(plc, 4, facet_line_ids, NULL);
//   }
//   {
//     viennagrid_element_id facet_line_ids[] = {lines[0], lines[2], lines[3], lines[4], lines[5], lines[6], lines[17], lines[18], lines[19]};
//     viennagrid_plc_facet_create(plc, 9, facet_line_ids, NULL);
//   }
//   {
//     viennagrid_element_id facet_line_ids[] = {lines[7], lines[8], lines[9], lines[10], lines[11]};
//     viennagrid_plc_facet_create(plc, 5, facet_line_ids, NULL);
//   }
//
//   {
//     viennagrid_element_id facet_line_ids[] = {lines[0], lines[2], lines[7], lines[12], lines[13], lines[24], lines[20], lines[25]};
//     viennagrid_plc_facet_create(plc, 8, facet_line_ids, NULL);
//   }
//   {
//     viennagrid_element_id facet_line_ids[] = {lines[3], lines[8], lines[13], lines[14]};
//     viennagrid_plc_facet_create(plc, 4, facet_line_ids, NULL);
//   }
//   {
//     viennagrid_element_id facet_line_ids[] = {lines[4], lines[9], lines[14], lines[15]};
//     viennagrid_plc_facet_create(plc, 4, facet_line_ids, NULL);
//   }
//   {
//     viennagrid_element_id facet_line_ids[] = {lines[5], lines[10], lines[15], lines[16]};
//     viennagrid_plc_facet_create(plc, 4, facet_line_ids, NULL);
//   }
//   {
//     viennagrid_element_id facet_line_ids[] = {lines[6], lines[11], lines[12], lines[16]};
//     viennagrid_plc_facet_create(plc, 4, facet_line_ids, NULL);
//   }
//
//   {
//     viennagrid_element_id facet_line_ids[] = {lines[17], lines[21], lines[24], lines[26]};
//     viennagrid_plc_facet_create(plc, 4, facet_line_ids, NULL);
//   }
//   {
//     viennagrid_element_id facet_line_ids[] = {lines[18], lines[22], lines[26], lines[27]};
//     viennagrid_plc_facet_create(plc, 4, facet_line_ids, NULL);
//   }
//   {
//     viennagrid_element_id facet_line_ids[] = {lines[19], lines[23], lines[25], lines[27]};
//     viennagrid_plc_facet_create(plc, 4, facet_line_ids, NULL);
//   }
//
//   return plc;
// }


