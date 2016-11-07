#include <string>
#include <sstream>

#include "common.hpp"
#include "templated_mesh.hpp"
#include "benchmark.hpp"


shared_ptr<TemplatedMesh> make_gummel_star_rot(BenchmarkConfig const & bench_cfg, bool coarse = false)
{
  MeshType gummel_geometry;
  {
    double gamma = (1.0+std::sqrt(5.0))/2.0;
    double gamma_2 = gamma*gamma;

    ElementType center = viennagrid::make_vertex( gummel_geometry, viennagrid::make_point(0, 0) );
    ElementType left = viennagrid::make_vertex( gummel_geometry, viennagrid::make_point(-std::sin(M_PI/5.0)/gamma_2, cos(M_PI/5.0)/gamma_2) );
    ElementType outer = viennagrid::make_vertex( gummel_geometry, viennagrid::make_point(0, 1) );
    ElementType right = viennagrid::make_vertex( gummel_geometry, viennagrid::make_point( std::sin(M_PI/5.0)/gamma_2, cos(M_PI/5.0)/gamma_2) );

    viennagrid::make_line(gummel_geometry, center, left);
    viennagrid::make_line(gummel_geometry, left, outer);
    viennagrid::make_line(gummel_geometry, outer, right);
    viennagrid::make_line(gummel_geometry, right, center);
  }

  MeshType gummel_mesh;
  gummel_mesh = refine_lines(gummel_geometry, bench_cfg.line_size());

  viennamesh::algorithm_handle gummel_mesher = context().make_algorithm("triangle_make_mesh");
  gummel_mesher.set_input("mesh", gummel_mesh);
  if (!coarse && (bench_cfg.cell_size() > 0))
  {
    gummel_mesher.set_input("cell_size", bench_cfg.cell_size());
  }
  if (!coarse && (bench_cfg.min_angle > 0))
  {
    gummel_mesher.set_input("min_angle", bench_cfg.min_angle);
  }

  gummel_mesher.set_input( "no_points_on_boundary", true );
  gummel_mesher.run();

  MeshType gummel_template_mesh = gummel_mesher.get_output<viennagrid_mesh>("mesh")();


  PointType min = viennagrid::make_point(-1, -1) * 1.01;
  PointType max = viennagrid::make_point(1, 1) * 1.01;

  shared_ptr<TemplatedMesh> templated_mesh(new TemplatedMesh(min, max));
  templated_mesh->add_instances_rotatation(5, gummel_template_mesh);
  templated_mesh->setup(bench_cfg.eps);

  return templated_mesh;
}

shared_ptr<TemplatedMesh> make_gummel_star_refl_rot(BenchmarkConfig const & bench_cfg, bool coarse = false)
{
  MeshType gummel_geometry;
  {
    double gamma = (1.0+std::sqrt(5.0))/2.0;
    double gamma_2 = gamma*gamma;

    ElementType center = viennagrid::make_vertex( gummel_geometry, viennagrid::make_point(0, 0) );
//     ElementType left = viennagrid::make_vertex( gummel_geometry, viennagrid::make_point(-std::sin(M_PI/5.0)/gamma_2, cos(M_PI/5.0)/gamma_2) );
    ElementType outer = viennagrid::make_vertex( gummel_geometry, viennagrid::make_point(0, 1) );
    ElementType right = viennagrid::make_vertex( gummel_geometry, viennagrid::make_point( std::sin(M_PI/5.0)/gamma_2, cos(M_PI/5.0)/gamma_2) );

    viennagrid::make_line(gummel_geometry, center, outer);
//     viennagrid::make_line(gummel_geometry, left, outer);
    viennagrid::make_line(gummel_geometry, outer, right);
    viennagrid::make_line(gummel_geometry, right, center);
  }

  MeshType gummel_mesh;
  gummel_mesh = refine_lines(gummel_geometry, bench_cfg.line_size());

  viennamesh::algorithm_handle gummel_mesher = context().make_algorithm("triangle_make_mesh");
  gummel_mesher.set_input("mesh", gummel_mesh);
  if (!coarse && (bench_cfg.cell_size() > 0))
  {
    gummel_mesher.set_input("cell_size", bench_cfg.cell_size());
  }
  if (!coarse && (bench_cfg.min_angle > 0))
  {
    gummel_mesher.set_input("min_angle", bench_cfg.min_angle);
  }

  gummel_mesher.set_input( "no_points_on_boundary", true );
  gummel_mesher.run();

  MeshType gummel_template_mesh = gummel_mesher.get_output<viennagrid_mesh>("mesh")();


  PointType min = viennagrid::make_point(-1, -1) * 1.01;
  PointType max = viennagrid::make_point(1, 1) * 1.01;

  shared_ptr<TemplatedMesh> templated_mesh(new TemplatedMesh(min, max));

  templated_mesh->add_template(gummel_template_mesh);

  Transformation mirror = Transformation::make_reflection_x(2);

  for (int i = 0; i != 5; ++i)
  {
    Transformation rotate = Transformation::make_rotation_2( i*2*M_PI/5 );
    templated_mesh->add_instance( 0, rotate, 0 );
    templated_mesh->add_instance( 0, composition(rotate, mirror), 0 );
  }

  templated_mesh->setup(bench_cfg.eps);

  return templated_mesh;
}





void write_dolfin_xml(MeshType const & mesh, std::string const & filename, double eps)
{
  std::ofstream mesh_file;
  mesh_file.open( (filename + ".xml").c_str() );

  std::map<ElementType, int> vertex_to_index_map;

  ElementRangeType vertices(mesh, 0);
  ElementRangeType triangles(mesh, 2);

  mesh_file << "<?xml version=\"1.0\"?>\n";
  mesh_file << "<dolfin xmlns:dolfin=\"http://fenicsproject.org\">\n";
  mesh_file << "  <mesh celltype=\"triangle\" dim=\"2\">\n";

  mesh_file << "    <vertices size=\"" << vertices.size() << "\">\n";
  int index = 0;
  for (ElementRangeIterator vit = vertices.begin(); vit != vertices.end(); ++vit, ++index)
  {
    PointType point = viennagrid::get_point(*vit);
    vertex_to_index_map[*vit] = index;

    mesh_file << "      <vertex index=\"" << index << "\" x=\"" << point[0] << "\" y=\"" << point[1] << "\" />\n";
  }
  mesh_file << "    </vertices>\n";

  mesh_file << "    <cells size=\"" << triangles.size() << "\">\n";
  index = 0;
  for (ElementRangeIterator tit = triangles.begin(); tit != triangles.end(); ++tit, ++index)
  {
    BoundaryElementRangeType bnd_vtx(*tit, 0);
    mesh_file << "      <triangle index=\"" << index << "\" v0=\"" << vertex_to_index_map[bnd_vtx[0]] <<
                                                    "\" v1=\"" << vertex_to_index_map[bnd_vtx[1]] <<
                                                    "\" v2=\"" << vertex_to_index_map[bnd_vtx[2]] << "\" />\n";
  }
  mesh_file << "    </cells>\n";

  mesh_file << "  </mesh>\n";
  mesh_file << "</dolfin>\n";





  double gamma = (1.0+std::sqrt(5.0))/2.0;
  double gamma_2 = gamma*gamma;

  PointType p[10];
  p[0] = viennagrid::make_point(0,1);
  p[1] = viennagrid::make_point(-std::sin(M_PI/5.0)/gamma_2, cos(M_PI/5.0)/gamma_2);

  for (int i = 1; i != 5; ++i)
  {
    double angle = -2.0*M_PI/5.0*i;
    Transformation rotation = Transformation::make_rotation_2(angle);

    p[2*i+0] = rotation(p[0]);
    p[2*i+1] = rotation(p[1]);
  }

  std::map<ElementType, int> vertex_to_boundary_map;

  for (ElementRangeIterator vit = vertices.begin(); vit != vertices.end(); ++vit)
  {
    if (!viennagrid::is_boundary(mesh, *vit))
      continue;

    PointType point = viennagrid::get_point(*vit);

    bool found = false;
    for (int i = 0; i != 10; ++i)
    {
      if (viennagrid::norm_2(point - p[i]) < eps)
      {
        vertex_to_boundary_map[*vit] = i;
        found = true;
        break;
      }
    }

    if (!found)
    {
      for (int i = 0; i != 10; ++i)
      {
        int j = (i+1)%10;

        viennagrid_bool is_inside;
        viennagrid_point_in_simplex_2(2, &point[0], &p[i][0], &p[j][0], eps, &is_inside);
        if (is_inside == VIENNAGRID_TRUE)
        {
          vertex_to_boundary_map[*vit] = 10+i;
          break;
        }
      }
    }
  }

  std::ofstream bnd_file;
  bnd_file.open( (filename + "_bnd.xml").c_str() );

  bnd_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
  bnd_file << "<dolfin xmlns:dolfin=\"http://fenicsproject.org\">\n";
  bnd_file << "  <mesh_function type=\"uint\" dim=\"0\" size=\"" << vertices.size() << "\">\n";

  for (ElementRangeIterator vit = vertices.begin(); vit != vertices.end(); ++vit)
  {
    std::map<ElementType, int>::const_iterator it = vertex_to_boundary_map.find(*vit);
    bnd_file << "    <entity index=\"" << vertex_to_index_map[*vit] << "\" value=\"" << (it != vertex_to_boundary_map.end() ? it->second : -1) << "\"/>\n";
  }

  bnd_file << "  </mesh_function>\n";
  bnd_file << "</dolfin>\n";

}



int main(int argc, char **argv)
{
  BenchmarkConfig bench_cfg(2);
  bench_cfg.from_args(argc, argv);

  shared_ptr<TemplatedMesh> rot_sym_mesh = make_gummel_star_rot(bench_cfg);
  shared_ptr<TemplatedMesh> refl_rot_sym_mesh = make_gummel_star_refl_rot(bench_cfg);
  shared_ptr<TemplatedMesh> coarse_mesh = make_gummel_star_refl_rot(bench_cfg, true);

  write_mesh(rot_sym_mesh->structure_instance, bench_cfg.output_filename + "_rotsym");
  write_dolfin_xml(rot_sym_mesh->structure_instance, bench_cfg.output_filename + "_rotsym", bench_cfg.eps);

  write_mesh(refl_rot_sym_mesh->structure_instance, bench_cfg.output_filename + "_reflrotsym");
  write_dolfin_xml(refl_rot_sym_mesh->structure_instance, bench_cfg.output_filename + "_reflrotsym", bench_cfg.eps);

  BenchmarkResult result;
  bench_conventional_2d(bench_cfg, result, coarse_mesh->structure_instance);

  write_mesh(result.mesh, bench_cfg.output_filename + "_nonsym");
  write_dolfin_xml(result.mesh, bench_cfg.output_filename + "_nonsym", bench_cfg.eps);

  return 0;
}










