#include "benchmark.hpp"
#include "3d_geometry.hpp"
#include <tclap/CmdLine.h>

bool BenchmarkConfig::from_args(int argc, char **argv)
{
  TCLAP::ValueArg<int> bench_count_input("b","bench-count", "Benchmark count", false, 1, "int");
  TCLAP::ValueArg<int> matrix_bench_count_input("m","matrix-bench-count", "Matrix benchmark count", false, -1, "int");
  TCLAP::ValueArg<double> eps_input("e","eps", "Benchmark count", false, 1e-6, "double");
  TCLAP::ValueArg<std::string> output_filename_input("o","output-file-prefix", "Output file prefix", false, "", "string");

  TCLAP::ValueArg<double> cell_size_input("c","cell-size", "Cell size", false, -1, "double");
  TCLAP::ValueArg<double> facet_size_input("f","facet-size", "Facet size", false, -1, "double");
  TCLAP::ValueArg<double> line_size_input("l","line-size", "Line size", false, -1, "double");

  TCLAP::ValueArg<double> min_angle_input("a","min-angle", "Min angle", false, -1, "double");
  TCLAP::ValueArg<double> min_angle_pi_input("","min-angle-pi-factor", "Min angle pi factor", false, -1, "double");
  TCLAP::ValueArg<double> radius_edge_ratio_input("r","radius-edge-ratio", "Radius edge radio", false, -1, "double");
  TCLAP::ValueArg<double> facet_min_angle_input("n","facet-min-angle", "Facet min angle", false, -1, "double");

  TCLAP::ValueArg<int> angle_histogram_bin_count_input("i","angle-histogram-bin-count", "Histrogram bin count", false, 18, "int");
  TCLAP::ValueArg<int> ratio_histogram_bin_count_input("j","ratio-histogram-bin-count", "Histrogram bin count", false, 20, "int");

  TCLAP::ValueArg<int> rotational_symmetry_order_input("y","rotational-symmetry-order", "Rotational symmetry order", false, -1, "int");
  TCLAP::ValueArg<std::string> seed_points_input("s","seed-points", "Seed Points", false, "", "string");
  TCLAP::UnlabeledValueArg<std::string> template_filename_input("template-filename", "File name of template mesh", false, "", "string");

  try
  {
    TCLAP::CmdLine cmd("Benchmark Reflection 2D", ' ', "1.0");

    cmd.add( bench_count_input );
    cmd.add( matrix_bench_count_input );
    cmd.add( eps_input );
    cmd.add( output_filename_input );

    cmd.add( cell_size_input );
    cmd.add( facet_size_input );
    cmd.add( line_size_input );

    cmd.add( min_angle_input );
    cmd.add( min_angle_pi_input );
    cmd.add( facet_min_angle_input );
    cmd.add( radius_edge_ratio_input );

    cmd.add( angle_histogram_bin_count_input );
    cmd.add( ratio_histogram_bin_count_input );

    cmd.add( rotational_symmetry_order_input );
    cmd.add( seed_points_input );
    cmd.add( template_filename_input );

    cmd.parse( argc, argv );


    bench_count = bench_count_input.getValue();
    matrix_bench_count = matrix_bench_count_input.getValue();

    eps = eps_input.getValue();
    output_filename = output_filename_input.getValue();

    cell_size_ = cell_size_input.getValue();
    facet_size_ = facet_size_input.getValue();
    line_size_ = line_size_input.getValue();

//     if ((cell_size > 0) && (facet_size < 0))
//     {
//       double tmp_ls = std::pow(12*cell_size/std::sqrt(2), 1.0/3.0);
//       facet_size = sqrt(3)/4*tmp_ls*tmp_ls * 0.5;
//     }
//
//     if ((cell_size > 0) && (line_size < 0))
//     {
//       line_size = std::pow(12*cell_size/std::sqrt(2), 1.0/3.0) * 0.5;
//     }


    min_angle = min_angle_input.getValue();
    if ((min_angle < 0) && (min_angle_pi_input.getValue() > 0))
      min_angle = M_PI / min_angle_pi_input.getValue();

    facet_min_angle = facet_min_angle_input.getValue();
    radius_edge_ratio = radius_edge_ratio_input.getValue();

    angle_hist_bin_count = angle_histogram_bin_count_input.getValue();
    ratio_hist_bin_count = ratio_histogram_bin_count_input.getValue();


    rotational_symmetry_order = rotational_symmetry_order_input.getValue();
    if (seed_points_input.getValue() != "")
    {
      std::vector<std::string> sps = split_string_brackets( seed_points_input.getValue(), ";" );
      for (std::size_t i = 0; i != sps.size(); ++i)
      {

        sps[i] = sps[i].substr(1, sps[i].size()-2);
        std::vector<std::string> sp = split_string_brackets( sps[i], ";" );

        if (sp.size() != 2)
        {
          std::cerr << "String to seed point conversion: an entry has no point and no region id: " << sps[i] << std::endl;
          return false;
        }

        PointType p = boost::lexical_cast<PointType>(sp[0]);
        viennagrid_region_id region_id = boost::lexical_cast<viennagrid_region_id>(sp[1]);

        seed_points.push_back( viennamesh::seed_point(p, region_id) );
      }
    }
    template_filename = template_filename_input.getValue();

    return true;
  }
  catch (TCLAP::ArgException &e)  // catch any exceptions
  {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return false;
  }
}


void BenchmarkResult::finalize(int angle_hist_bin_count, int ratio_hist_bin_count)
{
//   viennagrid_numeric best_ratio_2d = 1 / std::sqrt(3); // Delaunay Mesh Generation page 26
//   viennagrid_numeric best_ratio_3d = std::sqrt(6) / 4; // Delaunay Mesh Generation page 26
//   viennagrid_numeric best_ratio = (viennagrid::geometric_dimension(mesh) == 2) ? best_ratio_2d : best_ratio_3d;

  double angle_hist_start = 0;
  double angle_hist_end = M_PI/2;

  double ratio_hist_start = 1;
  double ratio_hist_end = 3;

  if (cell_count() > 0)
  {
    angle_field = min_angles(mesh);
    angle_statistic = statistics_from_quantity_field(mesh, angle_field, angle_hist_start, angle_hist_end, angle_hist_bin_count);
    angle_statistic.normalize();

    ratio_field = ::radius_edge_ratio(mesh);
    ratio_statistic = statistics_from_quantity_field(mesh, ratio_field, ratio_hist_start, ratio_hist_end, ratio_hist_bin_count);
    ratio_statistic.normalize();

    statistics_valid = true;
  }

  if (structure_instance_cell_count() > 0)
  {
    angle_field_templated = min_angles(templated_mesh->structure_instance);
    angle_statistic_templated = statistics_from_quantity_field(templated_mesh->structure_instance, angle_field_templated, angle_hist_start, angle_hist_end, angle_hist_bin_count);
    angle_statistic_templated.normalize();

    ratio_field_templated = ::radius_edge_ratio(templated_mesh->structure_instance);
    ratio_statistic_templated = statistics_from_quantity_field(templated_mesh->structure_instance, ratio_field_templated, ratio_hist_start, ratio_hist_end, ratio_hist_bin_count);
    ratio_statistic_templated.normalize();

    statistics_valid_templated = true;
  }


  SI_matrix = fromMesh(templated_mesh->structure_instance);
  templated_matrix = fromMesh(templated_mesh);
  mesh_matrix = fromMesh(mesh);
}


void BenchmarkResult::save_conventional(std::string const & filename) const
{
  viennagrid::io::vtk_writer<MeshType> writer;
  if (statistics_valid)
  {
    writer.add_scalar_data_on_cells( angle_field, "min_angle" );
    writer.add_scalar_data_on_cells( ratio_field, "radius_edge_ratio" );
  }
  writer( mesh, filename );
}

void BenchmarkResult::save_templated(std::string const & filename) const
{
  viennagrid::io::vtk_writer<MeshType> writer;
  if (statistics_valid_templated)
  {
    writer.add_scalar_data_on_cells( angle_field_templated, "min_angle" );
    writer.add_scalar_data_on_cells( ratio_field_templated, "radius_edge_ratio" );
    writer.add_scalar_data_on_cells( templated_mesh->template_index_field, "template_index" );
    writer.add_scalar_data_on_cells( templated_mesh->instance_index_field, "instance_index" );
  }
  writer( templated_mesh->structure_instance, filename );
}

void BenchmarkResult::save(std::string const & filename) const
{
  save_conventional(filename);
  save_templated(filename + "_templated");
}


void BenchmarkResult::print_conventional(BenchmarkConfig const & bench_cfg, std::ostream & stream) const
{
  stream << "cell count = " << cell_count() << "\n";
  stream << "time = " << time << "\n";
  stream << "size = " << size() << "\n";
  stream << "matrix size = " << mesh_matrix.memory_usage() << "\n";
  stream << "matrix time = " << time_mesh_matrix << "\n";

  stream << "Angle statistics\n";
  stream << angle_statistic << "\n";
  if (bench_cfg.min_angle > 0)
  {
    stream << "Under config value (" << bench_cfg.min_angle << ") = " << angle_statistic.histogram().sum_downwards(bench_cfg.min_angle) << "\n";
  }

  stream << "Radius-edge-ratio statistics\n";
  stream << ratio_statistic << "\n";
  if (bench_cfg.radius_edge_ratio > 0)
  {
    stream << "Over config value (" << bench_cfg.radius_edge_ratio << ") = " << ratio_statistic.histogram().sum_upwards(bench_cfg.radius_edge_ratio) << "\n";
  }
//   stream << "\n";
}

void BenchmarkResult::print_templated(BenchmarkConfig const & bench_cfg, std::ostream & stream) const
{
  stream << "Templated structure instance cell count = " << structure_instance_cell_count() << "\n";
  stream << "Templated time = " << time_templated() << "\n";
  stream << "Templated volume time = " << volume_time_templated << "\n";
  stream << "Templated facet time = " << facet_time_templated << "\n";
  stream << "Templated line time = " << line_time_templated << "\n";
  stream << "Templated size = " << size_templated() << "\n";
  stream << "Templated size SVB = " << size_templated_SVB() << "\n";
  stream << "Templated structure instance size = " << size_templated_structure_instance() << "\n";
  stream << "Templated structure instance matrix size = " << SI_matrix.memory_usage() << "\n";
  stream << "Templated structure instance matrix time = " << time_SI_matrix << "\n";
  stream << "Templated templated CSR matrix size = " << templated_matrix.memory_usage() << "\n";
  stream << "Templated templated CSR matrix time = " << time_templated_matrix << "\n";

  stream << "Templated angle statistics\n";
  stream << angle_statistic_templated << "\n";
  if (bench_cfg.min_angle > 0)
  {
    stream << "Under config value (" << bench_cfg.min_angle << ") = " << angle_statistic_templated.histogram().sum_downwards(bench_cfg.min_angle) << "\n";
  }

  stream << "Templated radius-edge-ratio statistics\n";
  stream << ratio_statistic_templated << "\n";
  if (bench_cfg.radius_edge_ratio > 0)
  {
    stream << "Over config value (" << bench_cfg.radius_edge_ratio << ") = " << ratio_statistic_templated.histogram().sum_upwards(bench_cfg.radius_edge_ratio) << "\n";
  }
//   stream << "\n";
}

void BenchmarkResult::print_benefits(std::ostream & stream) const
{
  stream << "benefit time              = " << benefit_time() << "\n";
  stream << "benefit time per cell     = " << benefit_time_per_cell() << "\n";

  stream << "benefit SI size              = " << structure_instance_benefit_size() << "\n";
  stream << "benefit SI size SVB          = " << structure_instance_benefit_size_SVB() << "\n";

  stream << "benefit remeshed size              = " << remeshed_benefit_size() << "\n";
  stream << "benefit remeshed size per cell     = " << remeshed_benefit_size_per_cell() << "\n";
  stream << "benefit remeshed size SVB          = " << remeshed_benefit_size_SVB() << "\n";
  stream << "benefit remeshed size SVB per cell = " << remeshed_benefit_size_SVB_per_cell() << "\n";

  stream << "CSR benefit SI size              = " << structure_instance_benefit_size_CSR() << "\n";
  stream << "CSR benefit SI size SVB          = " << structure_instance_benefit_size_SVB_CSR() << "\n";

  stream << "CSR benefit remeshed size              = " << remeshed_benefit_size_CSR() << "\n";
  stream << "CSR benefit remeshed size per cell     = " << remeshed_benefit_size_per_cell_CSR() << "\n";
  stream << "CSR benefit remeshed size SVB          = " << remeshed_benefit_size_SVB_CSR() << "\n";
  stream << "CSR benefit remeshed size SVB per cell = " << remeshed_benefit_size_SVB_per_cell_CSR() << "\n";

  stream << "templated CSR benefit SI size              = " << structure_instance_benefit_size_CSR_templated() << "\n";
  stream << "templated CSR benefit SI size SVB          = " << structure_instance_benefit_size_SVB_CSR_templated() << "\n";

  stream << "templated CSR benefit remeshed size              = " << remeshed_benefit_size_CSR_templated() << "\n";
  stream << "templated CSR benefit remeshed size per cell     = " << remeshed_benefit_size_per_cell_CSR_templated() << "\n";
  stream << "templated CSR benefit remeshed size SVB          = " << remeshed_benefit_size_SVB_CSR_templated() << "\n";
  stream << "templated CSR benefit remeshed size SVB per cell = " << remeshed_benefit_size_SVB_per_cell_CSR_templated() << "\n";

  stream << "matrix-vector product templated CSR vs remeshed CSR = " << time_templated_matrix/time_mesh_matrix << "\n";
  stream << "matrix-vector product templated CSR vs structure instance CSR = " << time_templated_matrix/time_SI_matrix << "\n";
}

void BenchmarkResult::print(BenchmarkConfig const & bench_cfg, std::ostream & stream) const
{
  print_conventional(bench_cfg, stream);
  stream << "\n";
  print_templated(bench_cfg, stream);
  stream << "\n";
  print_benefits(stream);
}



std::pair<PointType, PointType> oversized_bounding_box(MeshType const & mesh)
{
  std::pair<PointType, PointType> bb = viennagrid::bounding_box(mesh);
  viennagrid_numeric size = viennagrid::norm_2(bb.second - bb.first);
  PointType one( viennagrid::geometric_dimension(mesh), 1 );

  bb.first -= one * size / 100.0;
  bb.second += one * size / 100.0;
  return bb;
}



void bench_conventional_2d(BenchmarkConfig const & bench_cfg,
                           BenchmarkResult & bench_result,
                           MeshType const & structure_instance,
                           viennamesh::point_container const & hole_points)
{
  viennamesh::algorithm_handle boundary = context().make_algorithm("extract_boundary");
  boundary.set_input( "mesh", structure_instance );
  boundary.run();

  viennamesh::algorithm_handle geometry = context().make_algorithm("line_coarsening");
  geometry.set_default_source(boundary);
  geometry.set_input( "angle", 3.14 );
  geometry.run();

  viennamesh::algorithm_handle conventional_mesher = context().make_algorithm("triangle_make_mesh");
  conventional_mesher.link_input( "seed_points", boundary, "seed_points" );
  if (hole_points.empty())
    conventional_mesher.link_input( "hole_points", boundary, "hole_points" );
  else
  {
    for (std::size_t i = 0; i != hole_points.size(); ++i)
      conventional_mesher.push_back_input( "hole_points", hole_points[i] );
  }

  conventional_mesher.link_input( "mesh", geometry, "mesh" );
  if (bench_cfg.cell_size() > 0)
  {
    conventional_mesher.set_input("cell_size", bench_cfg.cell_size());
  }
  if (bench_cfg.min_angle > 0)
  {
    conventional_mesher.set_input("min_angle", bench_cfg.min_angle);
  }

  viennautils::Timer timer;
  timer.start();
  for (std::size_t k = 0; k != bench_cfg.bench_count; ++k)
  {
    conventional_mesher.run();
  }
  bench_result.time += timer.get() / bench_cfg.bench_count;
  bench_result.mesh = conventional_mesher.get_output<viennagrid_mesh>("mesh")();
}



shared_ptr<TemplatedMesh> make_templated_reflection_2d_x(MeshType template_mesh, viennagrid_numeric eps)
{
  std::pair<PointType, PointType> bb = oversized_bounding_box(template_mesh);
  shared_ptr<TemplatedMesh> templated_mesh(new TemplatedMesh( viennagrid::min(bb.first, -bb.second),
                                                             -viennagrid::min(bb.first, -bb.second)));

  templated_mesh->transformations.push_back( Transformation::make_identity(2) );
  templated_mesh->transformations.push_back( Transformation::make_reflection_x(2) );

  RegionRangeType regions(template_mesh);
  if (regions.empty())
  {
    templated_mesh->templates.push_back(template_mesh);

    MeshInstance mi0(templated_mesh);
    mi0.mesh_template_index = 0;
    mi0.transformation_index = 0;
    mi0.region_id = 0;
    templated_mesh->instances.push_back(mi0);

    MeshInstance mi1(templated_mesh);
    mi1.mesh_template_index = 0;
    mi1.transformation_index = 1;
    mi1.region_id = 0;
    templated_mesh->instances.push_back(mi1);
  }
  else
  {
    int index = 0;
    for (RegionRangeIterator rit = regions.begin(); rit != regions.end(); ++rit, ++index)
    {
      MeshType tmp;
      viennagrid::copy(*rit, tmp);
      templated_mesh->templates.push_back(tmp);

      MeshInstance mi0(templated_mesh);
      mi0.mesh_template_index = index;
      mi0.transformation_index = 0;
      mi0.region_id = (*rit).id();
      templated_mesh->instances.push_back(mi0);

      MeshInstance mi1(templated_mesh);
      mi1.mesh_template_index = index;
      mi1.transformation_index = 1;
      mi1.region_id = (*rit).id();
      templated_mesh->instances.push_back(mi1);
    }
  }

  templated_mesh->setup(eps);

  return templated_mesh;
}




MeshType bench_volumetric_mesh_generation_reflection_2d_x(BenchmarkConfig const & bench_cfg,
                                                          MeshType const & boundary_geometry,
                                                          double * time)
{
  viennamesh::algorithm_handle mesher = context().make_algorithm("triangle_make_mesh");

  mesher.set_input("mesh", boundary_geometry);
  if (time && (bench_cfg.cell_size() > 0))
  {
    mesher.set_input("cell_size", bench_cfg.cell_size());
  }
  if (time && (bench_cfg.min_angle > 0))
  {
    mesher.set_input("min_angle", bench_cfg.min_angle);
  }
  for (std::size_t i = 0; i != bench_cfg.seed_points.size(); ++i)
  {
    mesher.push_back_input( "seed_points", bench_cfg.seed_points[i] );
  }


  viennautils::Timer timer;
  timer.start();
  for (std::size_t k = 0; k != (time ? bench_cfg.bench_count : 1); ++k)
  {
    mesher.run();
  }
  double tmp = timer.get() / bench_cfg.bench_count;
  if (time) *time += tmp;
//   result.time_templated += timer.get() / bench_cfg.bench_count;

  return mesher.get_output<viennagrid_mesh>("mesh")();
}


BenchmarkResult bench_reflection_2d_x(BenchmarkConfig const & bench_cfg, MeshType const & boundary_geometry)

{
  BenchmarkResult result;

  MeshType template_mesh = bench_volumetric_mesh_generation_reflection_2d_x(bench_cfg, boundary_geometry, &result.volume_time_templated);
  result.templated_mesh = make_templated_reflection_2d_x(template_mesh, bench_cfg.eps);

  MeshType coarse_template_mesh = bench_volumetric_mesh_generation_reflection_2d_x(bench_cfg, boundary_geometry, NULL);
  boost::shared_ptr<TemplatedMesh> coarse_templated_mesh = make_templated_reflection_2d_x(coarse_template_mesh, bench_cfg.eps);

  bench_conventional_2d(bench_cfg, result, coarse_templated_mesh->structure_instance);

  return result;
}




shared_ptr<TemplatedMesh> make_templated_reflection_3d_xy(MeshType template_mesh, viennagrid_numeric eps)
{
  std::pair<PointType, PointType> bb = oversized_bounding_box(template_mesh);

  shared_ptr<TemplatedMesh> templated_mesh(new TemplatedMesh( viennagrid::min(bb.first, -bb.second),
                                                             -viennagrid::min(bb.first, -bb.second) ) );

  templated_mesh->transformations.push_back( Transformation::make_identity(3) );
  templated_mesh->transformations.push_back( Transformation::make_reflection_x(3) );
  templated_mesh->transformations.push_back( Transformation::make_reflection_y(3) );
  templated_mesh->transformations.push_back( Transformation::make_reflection_xy(3) );

  RegionRangeType regions(template_mesh);
  if (regions.empty())
  {
    templated_mesh->templates.push_back(template_mesh);

    MeshInstance mi0(templated_mesh);
    mi0.mesh_template_index = 0;
    mi0.transformation_index = 0;
    mi0.region_id = 0;
    templated_mesh->instances.push_back(mi0);

    MeshInstance mi1(templated_mesh);
    mi1.mesh_template_index = 0;
    mi1.transformation_index = 1;
    mi1.region_id = 0;
    templated_mesh->instances.push_back(mi1);

    MeshInstance mi2(templated_mesh);
    mi2.mesh_template_index = 0;
    mi2.transformation_index = 2;
    mi2.region_id = 0;
    templated_mesh->instances.push_back(mi2);

    MeshInstance mi3(templated_mesh);
    mi3.mesh_template_index = 0;
    mi3.transformation_index = 3;
    mi3.region_id = 0;
    templated_mesh->instances.push_back(mi3);
  }
  else
  {
    int index = 0;
    for (RegionRangeIterator rit = regions.begin(); rit != regions.end(); ++rit, ++index)
    {
      MeshType tmp;
      viennagrid::copy(*rit, tmp);
      templated_mesh->templates.push_back(tmp);

      MeshInstance mi0(templated_mesh);
      mi0.mesh_template_index = index;
      mi0.transformation_index = 0;
      mi0.region_id = (*rit).id();
      templated_mesh->instances.push_back(mi0);

      MeshInstance mi1(templated_mesh);
      mi1.mesh_template_index = index;
      mi1.transformation_index = 1;
      mi1.region_id = (*rit).id();
      templated_mesh->instances.push_back(mi1);

      MeshInstance mi2(templated_mesh);
      mi2.mesh_template_index = index;
      mi2.transformation_index = 2;
      mi2.region_id = (*rit).id();
      templated_mesh->instances.push_back(mi2);

      MeshInstance mi3(templated_mesh);
      mi3.mesh_template_index = index;
      mi3.transformation_index = 3;
      mi3.region_id = (*rit).id();
      templated_mesh->instances.push_back(mi3);
    }
  }

  templated_mesh->setup(eps);

  return templated_mesh;
}






void bench_conventional_3d(BenchmarkConfig const & bench_cfg,
                           BenchmarkResult & bench_result,
                           MeshType const & structure_instance)
{
  viennamesh::algorithm_handle boundary = context().make_algorithm("extract_boundary");
  boundary.set_input( "mesh", structure_instance );
  boundary.run();

  viennamesh::algorithm_handle geometry = context().make_algorithm("extract_plc_geometry");
  geometry.set_default_source(boundary);
  geometry.set_input( "coplanar_tolerance", 1e-8 );
  geometry.set_input( "colinear_tolerance", 1e-8 );
  geometry.run();

  viennamesh::algorithm_handle convert = context().make_algorithm("tetgen_convert");
  convert.set_default_source(geometry);
  convert.run();

  viennamesh::algorithm_handle conventional_mesher = context().make_algorithm("tetgen_make_mesh");
  conventional_mesher.link_input( "seed_points", boundary, "seed_points" );
  conventional_mesher.link_input( "hole_points", boundary, "hole_points" );
  conventional_mesher.link_input( "geometry", convert, "geometry" );
  if (bench_cfg.cell_size() > 0)
  {
    conventional_mesher.set_input("cell_size", bench_cfg.cell_size());
  }
  if (bench_cfg.min_angle > 0)
  {
    conventional_mesher.set_input("min_dihedral_angle", bench_cfg.min_angle);
  }
  if (bench_cfg.radius_edge_ratio > 0)
  {
    conventional_mesher.set_input("max_radius_edge_ratio", bench_cfg.radius_edge_ratio);
  }

  viennautils::Timer timer;
  timer.start();
  for (std::size_t k = 0; k != bench_cfg.bench_count; ++k)
  {
    conventional_mesher.run();
  }
  bench_result.time += timer.get() / bench_cfg.bench_count;
  bench_result.mesh = conventional_mesher.get_output<viennagrid_mesh>("mesh")();
}




BenchmarkResult bench_reflection_3d_xy(BenchmarkConfig const & bench_cfg, viennagrid_plc boundary_geometry)
{
  BenchmarkResult result;

  viennamesh::algorithm_handle convert = context().make_algorithm("tetgen_convert");
  convert.set_input("geometry", boundary_geometry);
  convert.run();


  viennamesh::algorithm_handle mesher = context().make_algorithm("tetgen_make_mesh");
  mesher.set_default_source(convert);
//   mesher.set_input("geometry", boundary_geometry);

  if (bench_cfg.cell_size() > 0)
  {
    mesher.set_input("cell_size", bench_cfg.cell_size());
  }
  if (bench_cfg.min_angle > 0)
  {
    mesher.set_input("min_dihedral_angle", bench_cfg.min_angle);
  }
  if (bench_cfg.radius_edge_ratio > 0)
  {
    mesher.set_input("max_radius_edge_ratio", bench_cfg.radius_edge_ratio);
  }
  for (std::size_t i = 0; i != bench_cfg.seed_points.size(); ++i)
  {
    mesher.push_back_input( "seed_points", bench_cfg.seed_points[i] );
  }

  viennautils::Timer timer;
  timer.start();
  for (std::size_t k = 0; k != bench_cfg.bench_count; ++k)
  {
    mesher.run();
  }
  result.volume_time_templated += timer.get() / bench_cfg.bench_count;
  result.templated_mesh = make_templated_reflection_3d_xy( mesher.get_output<viennagrid_mesh>("mesh")(), bench_cfg.eps );

  mesher.unset_input("cell_size");
  mesher.unset_input("min_dihedral_angle");
  mesher.unset_input("max_radius_edge_ratio");
  mesher.run();
  shared_ptr<TemplatedMesh> coarse_templated_mesh = make_templated_reflection_3d_xy( mesher.get_output<viennagrid_mesh>("mesh")(), bench_cfg.eps );

  bench_conventional_3d(bench_cfg, result, coarse_templated_mesh->structure_instance);

  return result;
}






shared_ptr<TemplatedMesh> make_templated_rotational_2d(MeshType template_mesh,
                                                       PointType const & min, PointType const & max,
                                                       int rotational_symmetry_order,
                                                       viennagrid_numeric eps)
{
  shared_ptr<TemplatedMesh> templated_mesh(new TemplatedMesh(min, max));

  double angle = 2*M_PI / rotational_symmetry_order;

  for (int i = 0; i != rotational_symmetry_order; ++i)
    templated_mesh->transformations.push_back( Transformation::make_rotation_2(angle*i) );
  templated_mesh->templates.push_back(template_mesh);

  for (int i = 0; i != rotational_symmetry_order; ++i)
  {
    MeshInstance mi0(templated_mesh);
    mi0.mesh_template_index = 0;
    mi0.transformation_index = i;
    mi0.region_id = 0;
    templated_mesh->instances.push_back(mi0);
  }

  templated_mesh->setup(eps);

  return templated_mesh;
}



MeshType refine_lines(MeshType const & mesh, viennagrid_numeric line_size)
{
  MeshType refined;

  ElementRangeType vertices(mesh, 0);
  std::vector<ElementType> new_vertices(vertices.size());

  int index = 0;
  for (ElementRangeIterator vit = vertices.begin(); vit != vertices.end(); ++vit, ++index)
    new_vertices[index] = viennagrid::make_vertex(refined, viennagrid::get_point(*vit));

  ElementRangeType lines(mesh, 1);
  for (ElementRangeIterator lit = lines.begin(); lit != lines.end(); ++lit)
  {
    ElementType vtx0 = viennagrid::vertices(*lit)[0];
    ElementType vtx1 = viennagrid::vertices(*lit)[1];

    PointType pt0 = viennagrid::get_point(vtx0);
    PointType pt1 = viennagrid::get_point(vtx1);

    PointType dir = pt1-pt0;
    int new_line_count = viennagrid::norm_2(dir)/line_size + 0.5;
    if (new_line_count < 1)
      new_line_count = 1;
    dir /= new_line_count;

    if (new_line_count > 1)
    {
      std::vector<ElementType> tmp_vtx(new_line_count-1);
      for (int i = 1; i != new_line_count; ++i)
        tmp_vtx[i-1] = viennagrid::make_vertex(refined, pt0+i*dir);

      viennagrid::make_line(refined, new_vertices[vtx0.id().index()], tmp_vtx[0]);
      for (int i = 0; i != new_line_count-2; ++i)
      {
        viennagrid::make_line(refined, tmp_vtx[i], tmp_vtx[i+1]);
      }
      viennagrid::make_line(refined, tmp_vtx[tmp_vtx.size()-1], new_vertices[vtx1.id().index()]);
    }
    else
    {
      viennagrid::make_line(refined, new_vertices[vtx0.id().index()], new_vertices[vtx1.id().index()]);
    }
  }

  return refined;
}







shared_ptr<TemplatedMesh> bench_rotational_2d_circle_impl(BenchmarkConfig const & bench_cfg, double radius,
                                                          double * line_time, double * volume_time)
{
  viennautils::Timer timer;

  PointType center = viennagrid::make_point(0,0);
  PointType top = viennagrid::make_point(0,radius);
  PointType left = rotate(top, -bench_cfg.angle()/2);
  PointType right = rotate(top, bench_cfg.angle()/2);


  MeshType boundary_geometry;
  ElementType center_vtx = viennagrid::make_vertex(boundary_geometry, center);
  ElementType left_vtx = viennagrid::make_vertex(boundary_geometry, left);
  ElementType right_vtx = viennagrid::make_vertex(boundary_geometry, right);

  viennagrid::make_line(boundary_geometry, center_vtx, right_vtx);
  viennagrid::make_line(boundary_geometry, right_vtx, left_vtx);
  viennagrid::make_line(boundary_geometry, left_vtx, center_vtx);

  MeshType boundary_mesh;

  if (line_time)
  {
    timer.start();
    for (std::size_t k = 0; k != bench_cfg.bench_count; ++k)
    {
      boundary_mesh = refine_lines(boundary_geometry, bench_cfg.line_size());
    }
    *line_time += timer.get() / bench_cfg.bench_count;
  }
  else
  {
    viennagrid::copy(boundary_geometry, boundary_mesh);
  }

  viennamesh::algorithm_handle mesher = context().make_algorithm("triangle_make_mesh");
  mesher.set_input("mesh", boundary_mesh);
  if (volume_time && (bench_cfg.cell_size() > 0))
  {
    mesher.set_input("cell_size", bench_cfg.cell_size());
  }
  if (volume_time && (bench_cfg.min_angle > 0))
  {
    mesher.set_input("min_angle", bench_cfg.min_angle);
  }

  mesher.set_input( "no_points_on_boundary", true );
  timer.start();
  for (std::size_t k = 0; k != (volume_time ? bench_cfg.bench_count : 1); ++k)
  {
    mesher.run();
  }
  if (volume_time)
    *volume_time += timer.get() / bench_cfg.bench_count;


  MeshType template_mesh = mesher.get_output<viennagrid_mesh>("mesh")();
  return make_templated_rotational_2d(template_mesh,
                                      -viennagrid::make_point(radius, radius) * 1.01,
                                      viennagrid::make_point(radius, radius) * 1.01,
                                      bench_cfg.rotational_symmetry_order, bench_cfg.eps);
}


BenchmarkResult bench_rotational_2d_circle(BenchmarkConfig const & bench_cfg, double radius)
{
  BenchmarkResult result;
  result.templated_mesh = bench_rotational_2d_circle_impl(bench_cfg, radius, &result.line_time_templated, &result.volume_time_templated);

  shared_ptr<TemplatedMesh> coarse_mesh = bench_rotational_2d_circle_impl(bench_cfg, radius, NULL, NULL);
  bench_conventional_2d(bench_cfg, result, coarse_mesh->structure_instance);

  return result;
}








shared_ptr<TemplatedMesh> make_templated_rotational_with_center(MeshType const & slice, MeshType const & inner,
                                                                BenchmarkConfig const & bench_cfg,
                                                                PointType const & min, PointType const & max)
{
  shared_ptr<TemplatedMesh> templated_mesh(new TemplatedMesh(min, max));

  for (int i = 0; i != bench_cfg.rotational_symmetry_order; ++i)
    templated_mesh->transformations.push_back( Transformation::make_rotation_2(bench_cfg.angle()*i) );

  templated_mesh->templates.push_back(inner);
  templated_mesh->templates.push_back(slice);

  {
    MeshInstance mi0(templated_mesh);
    mi0.mesh_template_index = 0;
    mi0.transformation_index = 0;
    mi0.region_id = 0;
    templated_mesh->instances.push_back(mi0);
  }

  for (int i = 0; i != bench_cfg.rotational_symmetry_order; ++i)
  {
    MeshInstance mi0(templated_mesh);
    mi0.mesh_template_index = 1;
    mi0.transformation_index = i;
    mi0.region_id = 0;
    templated_mesh->instances.push_back(mi0);
  }

  templated_mesh->setup(bench_cfg.eps);

  return templated_mesh;
}






shared_ptr<TemplatedMesh> bench_rotational_2d_circle_opt_impl(BenchmarkConfig const & bench_cfg, double radius, double inner_radius,
                                                              double * line_time, double * volume_time)
{
  viennautils::Timer timer;

  PointType center = viennagrid::make_point(0,0);
  PointType inner_top = viennagrid::make_point(0,inner_radius);
  PointType outer_top = viennagrid::make_point(0,radius);
  PointType inner_left = rotate(inner_top, -bench_cfg.angle()/2);
  PointType inner_right = rotate(inner_top, bench_cfg.angle()/2);
  PointType outer_left = rotate(outer_top, -bench_cfg.angle()/2);
  PointType outer_right = rotate(outer_top, bench_cfg.angle()/2);


  MeshType slice_geometry;
  {
    ElementType inner_left_vtx = viennagrid::make_vertex(slice_geometry, inner_left);
    ElementType inner_right_vtx = viennagrid::make_vertex(slice_geometry, inner_right);
    ElementType outer_left_vtx = viennagrid::make_vertex(slice_geometry, outer_left);
    ElementType outer_right_vtx = viennagrid::make_vertex(slice_geometry, outer_right);

    viennagrid::make_line(slice_geometry, inner_left_vtx, inner_right_vtx);
    viennagrid::make_line(slice_geometry, inner_right_vtx, outer_right_vtx);
    viennagrid::make_line(slice_geometry, outer_right_vtx, outer_left_vtx);
    viennagrid::make_line(slice_geometry, outer_left_vtx, inner_left_vtx);
  }

  MeshType inner_geometry;
  {
    std::vector<ElementType> vtx(bench_cfg.rotational_symmetry_order);
    for (int i = 0; i != bench_cfg.rotational_symmetry_order; ++i)
      vtx[i] = viennagrid::make_vertex(inner_geometry, rotate(inner_left, i*bench_cfg.angle()));

    for (int i = 0; i != bench_cfg.rotational_symmetry_order-1; ++i)
      viennagrid::make_line(inner_geometry, vtx[i], vtx[i+1]);
    viennagrid::make_line(inner_geometry, vtx[bench_cfg.rotational_symmetry_order-1], vtx[0]);
  }


  MeshType slice_boundary_mesh;
  MeshType inner_boundary_mesh;

  if (line_time)
  {
    timer.start();
    for (std::size_t k = 0; k != bench_cfg.bench_count; ++k)
    {
      slice_boundary_mesh = refine_lines(slice_geometry, bench_cfg.line_size());
      inner_boundary_mesh = refine_lines(inner_geometry, bench_cfg.line_size());
    }
    *line_time += timer.get() / bench_cfg.bench_count;
  }
  else
  {
    viennagrid::copy(slice_geometry, slice_boundary_mesh);
    viennagrid::copy(inner_geometry, inner_boundary_mesh);
  }


  viennamesh::algorithm_handle slice_mesher = context().make_algorithm("triangle_make_mesh");
  viennamesh::algorithm_handle inner_mesher = context().make_algorithm("triangle_make_mesh");
  slice_mesher.set_input("mesh", slice_boundary_mesh);
  inner_mesher.set_input("mesh", inner_boundary_mesh);
  if (volume_time && (bench_cfg.cell_size() > 0))
  {
    slice_mesher.set_input("cell_size", bench_cfg.cell_size());
    inner_mesher.set_input("cell_size", bench_cfg.cell_size());
  }
  if (volume_time && (bench_cfg.min_angle > 0))
  {
    slice_mesher.set_input("min_angle", bench_cfg.min_angle);
    inner_mesher.set_input("min_angle", bench_cfg.min_angle);
  }

  slice_mesher.set_input( "no_points_on_boundary", true );
  inner_mesher.set_input( "no_points_on_boundary", true );

  timer.start();
  for (std::size_t k = 0; k != (volume_time ? bench_cfg.bench_count : 1); ++k)
  {
    slice_mesher.run();
    inner_mesher.run();
  }
  if (volume_time)
    *volume_time += timer.get() / bench_cfg.bench_count;


  MeshType slice_template_mesh = slice_mesher.get_output<viennagrid_mesh>("mesh")();
  MeshType inner_template_mesh = inner_mesher.get_output<viennagrid_mesh>("mesh")();

  return make_templated_rotational_with_center(slice_template_mesh, inner_template_mesh,
                                               bench_cfg,
                                               -viennagrid::make_point(radius,radius) * 1.01,
                                               viennagrid::make_point(radius,radius) * 1.01);
}


BenchmarkResult bench_rotational_2d_circle_opt(BenchmarkConfig const & bench_cfg, double radius, double inner_radius)
{
  BenchmarkResult result;

  result.templated_mesh = bench_rotational_2d_circle_opt_impl(bench_cfg, radius, inner_radius, &result.line_time_templated, &result.volume_time_templated);

  shared_ptr<TemplatedMesh> coarse_mesh = bench_rotational_2d_circle_opt_impl(bench_cfg, radius, inner_radius, NULL, NULL);
  bench_conventional_2d(bench_cfg, result, coarse_mesh->structure_instance);

  return result;
}







shared_ptr<TemplatedMesh> make_bridge_2d(BridgeConfig const & bridge_cfg, BenchmarkConfig const & bench_cfg,
                                         double * line_time, double * volume_time,
                                         viennamesh::point_container * hole_points)
{
  viennautils::Timer timer;

  double beam_length = bridge_cfg.beam_length;
  double beam_thickness = bridge_cfg.beam_thickness;

  MeshType beam_geometry;
  {
    ElementType beam_left_bottom_vtx = viennagrid::make_vertex( beam_geometry, viennagrid::make_point(-beam_length/2, -beam_thickness/2) );
    ElementType beam_left_top_vtx = viennagrid::make_vertex( beam_geometry, viennagrid::make_point(-beam_length/2, beam_thickness/2) );
    ElementType beam_right_bottom_vtx = viennagrid::make_vertex( beam_geometry, viennagrid::make_point( beam_length/2, -beam_thickness/2) );
    ElementType beam_right_top_vtx = viennagrid::make_vertex( beam_geometry, viennagrid::make_point( beam_length/2, beam_thickness/2) );

    viennagrid::make_line(beam_geometry, beam_left_bottom_vtx, beam_left_top_vtx);
    viennagrid::make_line(beam_geometry, beam_left_top_vtx, beam_right_top_vtx);
    viennagrid::make_line(beam_geometry, beam_right_top_vtx, beam_right_bottom_vtx);
    viennagrid::make_line(beam_geometry, beam_right_bottom_vtx, beam_left_bottom_vtx);
  }

  MeshType beam_connector_geometry;
  {
    ElementType beam_left_bottom_vtx = viennagrid::make_vertex( beam_connector_geometry, viennagrid::make_point(-beam_thickness*sqrt(3)/2, -beam_thickness/2) );
    ElementType beam_left_top_vtx = viennagrid::make_vertex( beam_connector_geometry, viennagrid::make_point(-beam_thickness*sqrt(3)/2,  beam_thickness/2) );
    ElementType beam_middle_top_vtx = viennagrid::make_vertex( beam_connector_geometry, viennagrid::make_point( 0,  beam_thickness) );
    ElementType beam_right_bottom_vtx = viennagrid::make_vertex( beam_connector_geometry, viennagrid::make_point( beam_thickness*sqrt(3)/2, -beam_thickness/2) );
    ElementType beam_right_top_vtx = viennagrid::make_vertex( beam_connector_geometry, viennagrid::make_point( beam_thickness*sqrt(3)/2,  beam_thickness/2) );

    viennagrid::make_line(beam_connector_geometry, beam_left_bottom_vtx, beam_left_top_vtx);
    viennagrid::make_line(beam_connector_geometry, beam_left_top_vtx, beam_middle_top_vtx);
    viennagrid::make_line(beam_connector_geometry, beam_middle_top_vtx, beam_right_top_vtx);
    viennagrid::make_line(beam_connector_geometry, beam_right_top_vtx, beam_right_bottom_vtx);
    viennagrid::make_line(beam_connector_geometry, beam_right_bottom_vtx, beam_left_bottom_vtx);
  }



  MeshType beam_boundary_mesh;
  MeshType beam_connector_boundary_mesh;

  if (line_time)
  {
    timer.start();
    for (std::size_t k = 0; k != bench_cfg.bench_count; ++k)
    {
      beam_boundary_mesh = refine_lines(beam_geometry, bench_cfg.line_size());
      beam_connector_boundary_mesh = refine_lines(beam_connector_geometry, bench_cfg.line_size());
    }
    *line_time += timer.get() / bench_cfg.bench_count;
  }
  else
  {
    viennagrid::copy(beam_geometry, beam_boundary_mesh);
    viennagrid::copy(beam_connector_geometry, beam_connector_boundary_mesh);
  }


  viennamesh::algorithm_handle beam_mesher = context().make_algorithm("triangle_make_mesh");
  viennamesh::algorithm_handle beam_connector_mesher = context().make_algorithm("triangle_make_mesh");
  beam_mesher.set_input("mesh", beam_boundary_mesh);
  beam_connector_mesher.set_input("mesh", beam_connector_boundary_mesh);
  if (volume_time && (bench_cfg.cell_size() > 0))
  {
    beam_mesher.set_input("cell_size", bench_cfg.cell_size());
    beam_connector_mesher.set_input("cell_size", bench_cfg.cell_size());
  }
  if (volume_time && (bench_cfg.min_angle > 0))
  {
    beam_mesher.set_input("min_angle", bench_cfg.min_angle);
    beam_connector_mesher.set_input("min_angle", bench_cfg.min_angle);
  }

  beam_mesher.set_input( "no_points_on_boundary", true );
  beam_connector_mesher.set_input( "no_points_on_boundary", true );

  timer.start();
  for (std::size_t k = 0; k != (volume_time ? bench_cfg.bench_count : 1); ++k)
  {
    beam_mesher.run();
    beam_connector_mesher.run();
  }
  if (volume_time)
    *volume_time += timer.get() / bench_cfg.bench_count;


  MeshType beam_template_mesh = beam_mesher.get_output<viennagrid_mesh>("mesh")();
  MeshType beam_connector_template_mesh = beam_connector_mesher.get_output<viennagrid_mesh>("mesh")();




  double offset_connector_x = bridge_cfg.beam_length + bridge_cfg.beam_thickness * std::sqrt(3);

  PointType min = viennagrid::make_point( -bridge_cfg.beam_thickness*sqrt(3)/2, -bridge_cfg.beam_thickness ) * 1.01;
  PointType max = viennagrid::make_point( offset_connector_x*(bridge_cfg.instances+1) + bridge_cfg.beam_thickness*sqrt(3)/2, bridge_cfg.beam_length + bridge_cfg.beam_thickness) * 1.01;


  shared_ptr<TemplatedMesh> templated_mesh(new TemplatedMesh(min, max));

  Transformation rotate60  = Transformation::make_rotation_2( M_PI/3 );
  Transformation rotate120 = Transformation::make_rotation_2( 2*M_PI/3 );
  Transformation rotate180 = Transformation::make_rotation_2( M_PI );

  templated_mesh->templates.push_back(beam_template_mesh);
  templated_mesh->templates.push_back(beam_connector_template_mesh);




  templated_mesh->add_instance(1, Transformation::make_translation(viennagrid::make_point(0, 0)), 0 );
  templated_mesh->add_instance(1, Transformation::make_translation(viennagrid::make_point(offset_connector_x, 0)), 0 );
  templated_mesh->add_instance(1, composition(Transformation::make_translation(viennagrid::make_point(offset_connector_x/2, std::sin(M_PI/3)*offset_connector_x)), rotate180), 0 );

  templated_mesh->add_instance(0, Transformation::make_translation(viennagrid::make_point(offset_connector_x/2, 0)), 0 );
  templated_mesh->add_instance(0, composition(Transformation::make_translation(viennagrid::make_point(std::cos(M_PI/3)*offset_connector_x/2+offset_connector_x/2, std::sin(M_PI/3)*offset_connector_x/2)), rotate60), 0 );
  templated_mesh->add_instance(0, composition(Transformation::make_translation(viennagrid::make_point(std::cos(2*M_PI/3)*offset_connector_x/2+offset_connector_x/2, std::sin(2*M_PI/3)*offset_connector_x/2)), rotate120), 0 );

  if (hole_points)
    hole_points->push_back( viennagrid::make_point(offset_connector_x/2, std::sin(M_PI/3)*offset_connector_x/2) );

  for (int i = 0; i != bridge_cfg.instances; ++i)
  {
    double offset_x = offset_connector_x*(i+1);

    templated_mesh->add_instance(1, Transformation::make_translation(viennagrid::make_point(offset_x+offset_connector_x, 0)), 0 );
    templated_mesh->add_instance(1, composition(Transformation::make_translation(viennagrid::make_point(offset_connector_x/2+offset_x, std::sin(M_PI/3)*offset_connector_x)), rotate180), 0 );

    templated_mesh->add_instance(0, Transformation::make_translation(viennagrid::make_point(offset_x+offset_connector_x/2, 0)), 0 );
    templated_mesh->add_instance(0, Transformation::make_translation(viennagrid::make_point(offset_x, std::sin(M_PI/3)*offset_connector_x)), 0 );
    templated_mesh->add_instance(0, composition(Transformation::make_translation(viennagrid::make_point(offset_x+std::cos(M_PI/3)*offset_connector_x/2+offset_connector_x/2, std::sin(M_PI/3)*offset_connector_x/2)), rotate60), 0 );
    templated_mesh->add_instance(0, composition(Transformation::make_translation(viennagrid::make_point(offset_x+std::cos(2*M_PI/3)*offset_connector_x/2+offset_connector_x/2, std::sin(2*M_PI/3)*offset_connector_x/2)), rotate120), 0 );

    if (hole_points)
    {
      hole_points->push_back( viennagrid::make_point(offset_x+offset_connector_x/2, std::sin(M_PI/3)*offset_connector_x/2) );
      hole_points->push_back( viennagrid::make_point(offset_x, std::sin(M_PI/3)*offset_connector_x/2) );
    }
  }


  templated_mesh->setup(bench_cfg.eps);

  return templated_mesh;
}


BenchmarkResult bench_bridge_2d(BridgeConfig const & bridge_cfg, BenchmarkConfig const & bench_cfg)
{
  BenchmarkResult result;

  result.templated_mesh = make_bridge_2d(bridge_cfg, bench_cfg, &result.line_time_templated, &result.volume_time_templated, NULL);

  viennamesh::point_container hole_points;
  shared_ptr<TemplatedMesh> coarse_mesh = make_bridge_2d(bridge_cfg, bench_cfg, NULL, NULL, &hole_points);
  bench_conventional_2d(bench_cfg, result, coarse_mesh->structure_instance, hole_points);

  return result;
}










TemplatedSurface make_3d_surface_mesh(viennagrid_plc plc, viennagrid_numeric eps)
{
  TemplatedSurface surface;

  viennagrid_dimension geo_dim;
  viennagrid_plc_geometric_dimension_get(plc, &geo_dim);

  viennagrid_element_id facets_begin;
  viennagrid_element_id facets_end;
  viennagrid_plc_elements_get(plc, 2, &facets_begin, &facets_end);

  for (viennagrid_element_id fid = facets_begin; fid != facets_end; ++fid)
  {
    viennagrid_element_id * facet_vertices_begin;
    viennagrid_element_id * facet_vertices_end;
    viennagrid_plc_boundary_elements(plc, fid, 0, &facet_vertices_begin, &facet_vertices_end);
    viennagrid_int facet_vertex_count = facet_vertices_end-facet_vertices_begin;

    viennagrid_element_id * facet_lines_begin;
    viennagrid_element_id * facet_lines_end;
    viennagrid_plc_boundary_elements(plc, fid, 1, &facet_lines_begin, &facet_lines_end);
    viennagrid_int facet_line_count = facet_lines_end-facet_lines_begin;

    std::pair<Transformation, Transformation> T = Transformation::make_facet_projection(plc, fid, eps);


    TemplatedFacet facet;
    facet.from_2d = T.second;

    facet.facet_template->geometry.init_points(facet_vertex_count);
    facet.facet_template->geometry.init_segments(facet_line_count);


    std::map<viennagrid_element_id, int> global_to_local_mapping;

    int index = 0;
    for (viennagrid_element_id * vit = facet_vertices_begin; vit != facet_vertices_end; ++vit, ++index)
    {
      PointType p3d = get_point(plc, *vit);
      PointType p2d = T.first(p3d);
      facet.facet_template->geometry.pointlist[2*index+0] = p2d[0];
      facet.facet_template->geometry.pointlist[2*index+1] = p2d[1];

      global_to_local_mapping[*vit] = index;
    }

    index = 0;
    for (viennagrid_element_id * lit = facet_lines_begin; lit != facet_lines_end; ++lit, ++index)
    {
      viennagrid_element_id * line_vertices_begin;
      viennagrid_element_id * line_vertices_end;
      viennagrid_plc_boundary_elements(plc, *lit, 0, &line_vertices_begin, &line_vertices_end);

      facet.facet_template->geometry.segmentlist[2*index+0] = global_to_local_mapping[ line_vertices_begin[0] ];
      facet.facet_template->geometry.segmentlist[2*index+1] = global_to_local_mapping[ line_vertices_begin[1] ];
    }



    viennagrid_int facet_hole_points_count;
    viennagrid_numeric * facet_hole_point_coords;
    viennagrid_plc_facet_hole_points_get(plc, fid, &facet_hole_points_count, &facet_hole_point_coords);

    facet.facet_template->geometry.init_hole_points(facet_hole_points_count);
    for (int i = 0; i != facet_hole_points_count; ++i)
    {
      PointType hp3d(geo_dim, facet_hole_point_coords+geo_dim*i);
      PointType hp2d = T.first(hp3d);

      facet.facet_template->geometry.holelist[2*i+0] = hp2d[0];
      facet.facet_template->geometry.holelist[2*i+1] = hp2d[1];
    }


    surface.facets.push_back( facet );
  }

  return surface;
}


























void bench_line_refine(BenchmarkConfig const & bench_cfg, viennagrid_plc & plc, double * time)
{
  viennautils::Timer timer;

  if (time && (bench_cfg.line_size() > 0))
  {
    viennagrid_plc refined_plc;
    viennagrid_plc_create(&refined_plc);

    timer.start();
    for (std::size_t k = 0; k != bench_cfg.bench_count; ++k)
    {
      viennagrid_plc_line_refine(plc, refined_plc, bench_cfg.line_size());
    }
    double tmp = timer.get() / bench_cfg.bench_count;

//     std::cout << "  bench line: " << tmp << std::endl;
    if (time) *time += tmp;

    viennagrid_plc_release(plc);
    plc = refined_plc;
  }
}

MeshType bench_surface_generation(BenchmarkConfig const & bench_cfg, TemplatedSurface & surface,
                                  PointType const & min, PointType const & max,
                                  viennagrid_numeric * time)
{
  viennautils::Timer timer;
  std::ostringstream surface_meshing_options;
  surface_meshing_options << "zpQYY";

  if (time && (bench_cfg.facet_min_angle > 0))
    surface_meshing_options << "q" << bench_cfg.facet_min_angle / M_PI * 180.0;

  if (time && (bench_cfg.facet_size() > 0))
    surface_meshing_options << "a" << std::fixed << std::showpoint << bench_cfg.facet_size();

  timer.start();
  for (std::size_t k = 0; k != (time ? bench_cfg.bench_count : 1); ++k)
  {
    surface.reset_meshes();
    surface.generate_meshes( surface_meshing_options.str() );
  }
  double tmp = timer.get() / bench_cfg.bench_count;

//   std::cout << "  bench surface: " << tmp << std::endl;
  if (time) *time += tmp;

  return surface.compose(min, max, bench_cfg.eps);
}



MeshType bench_volume_generation(BenchmarkConfig const & bench_cfg, viennagrid_plc geometry, MeshType const & surface_mesh, viennagrid_numeric * time)
{
  viennautils::Timer timer;

  viennamesh::algorithm_handle mesher = context().make_algorithm("tetgen_make_mesh");
  mesher.set_input( "geometry", surface_mesh );

  viennagrid_int seed_point_count;
  viennagrid_numeric * seed_point_coords;
  viennagrid_int * seed_point_regions;
  viennagrid_plc_seed_points_get(geometry, &seed_point_count, &seed_point_coords, &seed_point_regions);

  for (viennagrid_int i = 0; i != seed_point_count; ++i)
  {
    mesher.push_back_input("seed_points", viennamesh::seed_point(PointType(3, seed_point_coords + 3*i), seed_point_regions[i]));
  }

  mesher.set_input("forbid_steiner_points_on_faces", true);

  if (time && (bench_cfg.cell_size() > 0))
  {
    mesher.set_input("cell_size", bench_cfg.cell_size());
  }
  if (time && (bench_cfg.min_angle > 0))
  {
    mesher.set_input("min_dihedral_angle", bench_cfg.min_angle);
  }
  if (time && (bench_cfg.radius_edge_ratio > 0))
  {
    mesher.set_input("max_radius_edge_ratio", bench_cfg.radius_edge_ratio);
  }


  timer.start();
  for (std::size_t k = 0; k != (time ? bench_cfg.bench_count : 1); ++k)
  {
    mesher.run();
  }
  double tmp = timer.get() / bench_cfg.bench_count;

//   std::cout << "  bench volume: " << tmp << std::endl;
  if (time) *time += tmp;

  return mesher.get_output<viennagrid_mesh>("mesh")();
}




shared_ptr<TemplatedMesh> make_rotational_3d_tsv(BenchmarkConfig const & bench_cfg, TSVconfig const & tsv_cfg,
                                                 viennagrid_numeric * line_time, viennagrid_numeric * facet_time, viennagrid_numeric * volume_time)
{
  viennagrid_plc plc = make_rotational_3d_geometry(tsv_cfg, bench_cfg.eps);
  bench_line_refine(bench_cfg, plc, line_time);

  TemplatedSurface surface = make_3d_surface_mesh(plc, bench_cfg.eps);


  Transformation rot = find_transform(surface.facets[0].facet_template->geometry,
                                      surface.facets[1].facet_template->geometry, bench_cfg.eps);
  surface.facets[1].facet_template = surface.facets[0].facet_template;
  surface.facets[1].from_2d = composition(surface.facets[1].from_2d, rot);


  MeshType surface_mesh = bench_surface_generation(bench_cfg, surface, tsv_cfg.min(), tsv_cfg.max(), facet_time);

  MeshType template_mesh = bench_volume_generation(bench_cfg, plc, surface_mesh, volume_time);


  viennagrid_plc_release(plc);

  shared_ptr<TemplatedMesh> templated_mesh( new TemplatedMesh(tsv_cfg.min(), tsv_cfg.max()) );
  templated_mesh->add_instances_rotatation_z(tsv_cfg.rotational_symmetry_order, template_mesh);
  templated_mesh->setup(bench_cfg.eps);
  return templated_mesh;
}

BenchmarkResult bench_rotational_3d_tsv(BenchmarkConfig const & bench_cfg, TSVconfig const & cfg)
{
  BenchmarkResult result;

  result.templated_mesh = make_rotational_3d_tsv(bench_cfg, cfg, &result.line_time_templated, &result.facet_time_templated, &result.volume_time_templated);
  shared_ptr<TemplatedMesh> coarse_templated_mesh = make_rotational_3d_tsv(bench_cfg, cfg, NULL, NULL, NULL);
  bench_conventional_3d(bench_cfg, result, coarse_templated_mesh->structure_instance);

  return result;
}






shared_ptr<TemplatedMesh> make_rotational_3d_tsv_nonsym(BenchmarkConfig const & bench_cfg, TSVconfig const & tsv_cfg,
                                                        viennagrid_numeric * line_time, viennagrid_numeric * facet_time, viennagrid_numeric * volume_time)
{
  viennagrid_plc sym_plc = make_rotational_3d_geometry_nonsym(tsv_cfg, bench_cfg.eps);
  viennagrid_plc nonsym_plc = make_rotational_3d_geometry_nonsym_middle(tsv_cfg);


  bench_line_refine(bench_cfg, sym_plc, line_time);
  bench_line_refine(bench_cfg, nonsym_plc, line_time);


  TemplatedSurface sym_surface = make_3d_surface_mesh(sym_plc, bench_cfg.eps);
  TemplatedSurface nonsym_surface = make_3d_surface_mesh(nonsym_plc, bench_cfg.eps);

  int sym_inner_id = -1;
  {
    for (std::size_t i = 0; i != sym_surface.facets.size(); ++i)
    {
      Transformation rot = find_transform(sym_surface.facets[i].facet_template->geometry,
                                          nonsym_surface.facets[0].facet_template->geometry, bench_cfg.eps);

      if (rot.valid())
      {
        sym_inner_id = i;
        break;
      }
    }
  }


  TemplatedFacet & shared_facet_sym = sym_surface.facets[sym_inner_id];
  for (int i = 0; i != tsv_cfg.rotational_symmetry_order; ++i)
  {
    Transformation rot = find_transform(shared_facet_sym.facet_template->geometry,
                                        nonsym_surface.facets[i].facet_template->geometry, bench_cfg.eps);

    nonsym_surface.facets[i].facet_template = shared_facet_sym.facet_template;
    nonsym_surface.facets[i].from_2d = composition(nonsym_surface.facets[i].from_2d, rot);
  }

  MeshType sym_surface_mesh = bench_surface_generation(bench_cfg, sym_surface, tsv_cfg.min(), tsv_cfg.max(), facet_time);
  MeshType nonsym_surface_mesh = bench_surface_generation(bench_cfg, nonsym_surface, tsv_cfg.min(), tsv_cfg.max(), facet_time);


  MeshType sym_template_mesh = bench_volume_generation(bench_cfg, sym_plc, sym_surface_mesh, volume_time);
  MeshType nonsym_template_mesh = bench_volume_generation(bench_cfg, nonsym_plc, nonsym_surface_mesh, volume_time);


  viennagrid_plc_release(sym_plc);
  viennagrid_plc_release(nonsym_plc);



  shared_ptr<TemplatedMesh> templated_mesh( new TemplatedMesh(tsv_cfg.min(), tsv_cfg.max()) );
  templated_mesh->add_instances_rotatation_z(tsv_cfg.rotational_symmetry_order, sym_template_mesh);
  templated_mesh->add_trivial_instance(nonsym_template_mesh, 0);
  templated_mesh->setup(bench_cfg.eps);

  return templated_mesh;
}

BenchmarkResult bench_rotational_3d_tsv_nonsym(BenchmarkConfig const & bench_cfg, TSVconfig const & tsv_cfg)
{
  BenchmarkResult result;

  result.templated_mesh = make_rotational_3d_tsv_nonsym(bench_cfg, tsv_cfg, &result.line_time_templated, &result.facet_time_templated, &result.volume_time_templated);
  shared_ptr<TemplatedMesh> coarse_templated_mesh = make_rotational_3d_tsv(bench_cfg, tsv_cfg, NULL, NULL, NULL);
  bench_conventional_3d(bench_cfg, result, coarse_templated_mesh->structure_instance);

  return result;
}











shared_ptr<TemplatedMesh> make_multi_tsv(BenchmarkConfig const & bench_cfg, TSVconfig const & tsv_cfg,
                                         viennagrid_numeric * line_time, viennagrid_numeric * facet_time, viennagrid_numeric * volume_time)
{
  std::vector<PointType> centers;

  centers.push_back( viennagrid::make_point(0,0,0) );
  centers.push_back( viennagrid::make_point(0,1,0) * tsv_cfg.distance );
  centers.push_back( viennagrid::make_point(std::cos(M_PI/6),std::sin(M_PI/6),0) * tsv_cfg.distance );
  centers.push_back( viennagrid::make_point(std::cos(M_PI/6),-std::sin(M_PI/6),0) * tsv_cfg.distance );
  centers.push_back( viennagrid::make_point(0,-1,0) * tsv_cfg.distance );
  centers.push_back( viennagrid::make_point(-std::cos(M_PI/6),-std::sin(M_PI/6),0) * tsv_cfg.distance );
  centers.push_back( viennagrid::make_point(-std::cos(M_PI/6),std::sin(M_PI/6),0) * tsv_cfg.distance );



  viennagrid_plc tsv_slice_plc = make_rotational_3d_geometry(tsv_cfg, bench_cfg.eps);
  viennagrid_plc trunk_plc = make_multi_tsv_geometry(tsv_cfg, centers, bench_cfg.eps);


  bench_line_refine(bench_cfg, tsv_slice_plc, line_time);
  bench_line_refine(bench_cfg, trunk_plc, line_time);


  TemplatedSurface tsv_slice_surface = make_3d_surface_mesh(tsv_slice_plc, bench_cfg.eps);
  TemplatedSurface trunk_surface = make_3d_surface_mesh(trunk_plc, bench_cfg.eps);


  Transformation rot = find_transform(tsv_slice_surface.facets[0].facet_template->geometry,
                                      tsv_slice_surface.facets[1].facet_template->geometry, bench_cfg.eps);
  tsv_slice_surface.facets[1].facet_template = tsv_slice_surface.facets[0].facet_template;
  tsv_slice_surface.facets[1].from_2d = composition(tsv_slice_surface.facets[1].from_2d, rot);


  int interface_facet_id = -1;
  {
    for (std::size_t i = 0; i != tsv_slice_surface.facets.size(); ++i)
    {
      Transformation rot = find_transform(tsv_slice_surface.facets[i].facet_template->geometry,
                                          trunk_surface.facets[0].facet_template->geometry, bench_cfg.eps);

      if (rot.valid())
      {
        interface_facet_id = i;
        break;
      }
    }
  }


  TemplatedFacet & shared_facet_sym = tsv_slice_surface.facets[interface_facet_id];
  for (int i = 0; i != centers.size()*tsv_cfg.rotational_symmetry_order; ++i)
  {
    bool swap = i % tsv_cfg.rotational_symmetry_order == tsv_cfg.rotational_symmetry_order-1; // swapping last rotational facet to ensure same facet orientation in the structure instance

    Transformation rot = find_transform(shared_facet_sym.facet_template->geometry,
                                        trunk_surface.facets[i].facet_template->geometry, bench_cfg.eps, swap);

    trunk_surface.facets[i].facet_template = shared_facet_sym.facet_template;
    trunk_surface.facets[i].from_2d = composition(trunk_surface.facets[i].from_2d, rot);
  }

  MeshType tsv_slice_surface_mesh = bench_surface_generation(bench_cfg, tsv_slice_surface, tsv_cfg.multi_min(), tsv_cfg.multi_max(), facet_time);
  MeshType trunk_surface_mesh = bench_surface_generation(bench_cfg, trunk_surface, tsv_cfg.multi_min(), tsv_cfg.multi_max(), facet_time);


  MeshType tsv_slice_mesh = bench_volume_generation(bench_cfg, tsv_slice_plc, tsv_slice_surface_mesh, volume_time);
  MeshType trunk_mesh = bench_volume_generation(bench_cfg, trunk_plc, trunk_surface_mesh, volume_time);


  viennagrid_plc_release(tsv_slice_plc);
  viennagrid_plc_release(trunk_plc);


  shared_ptr<TemplatedMesh> templated_mesh( new TemplatedMesh(tsv_cfg.multi_min(), tsv_cfg.multi_max()) );
  templated_mesh->add_instances_rotatation_z(tsv_cfg.rotational_symmetry_order, tsv_slice_mesh, centers);
  templated_mesh->add_trivial_instance(trunk_mesh, 0);
  templated_mesh->setup(bench_cfg.eps);

  return templated_mesh;
}


BenchmarkResult bench_multi_tsv(BenchmarkConfig const & bench_cfg, TSVconfig const & tsv_cfg)
{
  BenchmarkResult result;

  result.templated_mesh = make_multi_tsv(bench_cfg, tsv_cfg, &result.line_time_templated, &result.facet_time_templated, &result.volume_time_templated);
  shared_ptr<TemplatedMesh> coarse_templated_mesh = make_multi_tsv(bench_cfg, tsv_cfg, NULL, NULL, NULL);
  bench_conventional_3d(bench_cfg, result, coarse_templated_mesh->structure_instance);

  return result;
}

