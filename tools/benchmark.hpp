#ifndef DISS_BENCHMARK_HPP_
#define DISS_BENCHMARK_HPP_

#include "common.hpp"
#include "templated_mesh.hpp"
#include "matrix.hpp"


struct BenchmarkConfig
{
  BenchmarkConfig(int dim_) : dim(dim_), bench_count(1), matrix_bench_count(-1), line_size_(-1), facet_size_(-1), cell_size_(-1), min_angle(-1), radius_edge_ratio(-1), facet_min_angle(-1), eps(1e-6) {}

  int dim;

  int bench_count;
  int matrix_bench_count;

  double eps;
  std::string output_filename;

  double line_size_;
  double facet_size_;
  double cell_size_;


  double cell_size() const
  {
    return cell_size_;
  }

  double facet_size() const
  {
    if ((cell_size_ > 0) && (facet_size_ < 0))
    {
      double tmp_ls = std::pow(12*cell_size_/std::sqrt(2), 1.0/3.0);
      return sqrt(3)/4*tmp_ls*tmp_ls * 0.5;
    }
    else
      return facet_size_;
  }

  double line_size() const
  {
    if (dim == 2)
    {
      return sqrt(cell_size() * 4.0 / sqrt(3.0)) * 0.75;
    }
    else
    {
      if ((cell_size_ > 0) && (line_size_ < 0))
      {
        return std::pow(12*cell_size_/std::sqrt(2), 1.0/3.0) * 0.5;
      }
      else
        return line_size_;
    }
  }

  double inner_radius(double radius = 1) const
  {
    return std::min(line_size()*rotational_symmetry_order / (2*M_PI), radius * 0.5);
  }

  double min_angle;
  double radius_edge_ratio;
  double facet_min_angle;

  int angle_hist_bin_count;
  int ratio_hist_bin_count;

  std::string template_filename;
  int rotational_symmetry_order;
  viennamesh::seed_point_container seed_points;

  double angle() const
  {
    return 2*M_PI / rotational_symmetry_order;
  }

  bool output_active() const
  {
    return !output_filename.empty();
  }

  std::string output_postfix() const
  {
    std::stringstream ss;

    if (cell_size_ > 0)
      ss << "_c" << cell_size_;
    if (facet_size_ > 0)
      ss << "_f" << facet_size_;
    if (line_size_ > 0)
      ss << "_l" << line_size_;

    if (min_angle > 0)
      ss << "_a" << min_angle;
    if (radius_edge_ratio > 0)
      ss << "_r" << radius_edge_ratio;
    if (facet_min_angle > 0)
      ss << "_n" << facet_min_angle;

    if (rotational_symmetry_order > 0)
      ss << "_y" << rotational_symmetry_order;

    return ss.str();
  }

  std::string full_output_filename() const
  {
    return output_filename + output_postfix();
  }

  bool from_args(int argc, char **argv);

  void print(std::ostream & stream = std::cout) const
  {
    stream << "benchmark count = " << bench_count << "\n";
    stream << "line size = " << line_size() << "   " << line_size_ << "\n";
    stream << "facet size = " << facet_size() << "   " << facet_size_ << "\n";
    stream << "cell size = " << cell_size() << "\n";
    stream << "min angle = " << min_angle << "\n";
    stream << "radius edge ratio = " << radius_edge_ratio << "\n";
    stream << "facet min angle = " << facet_min_angle << "\n";
    stream << "rotational symmetry order = " << rotational_symmetry_order << "\n";
    stream << "inner radius = " << line_size()*rotational_symmetry_order / (2*M_PI) << "\n";
    stream << "\n";
  }
};


struct TSVconfig
{
  double radius;
  double total_height;

  double hole_radius;
  double hole_height;

  double layer_thickness;

  double nonsym_inner_radius;

  double trunk_height() const { return total_height - hole_height - layer_thickness; }
  double body_height() const { return hole_height - 2*layer_thickness; }

  int rotational_symmetry_order;
  double angle() const { return 2*M_PI / rotational_symmetry_order; }

  PointType min() const { return viennagrid::make_point(-radius, -radius, -1) * 1.01; }
  PointType max() const { return viennagrid::make_point(radius, radius, total_height) * 1.01; }

  double distance;
  double trunk_size;

  PointType multi_min() const { return viennagrid::make_point(-trunk_size, -trunk_size, -1) * 1.01; }
  PointType multi_max() const { return viennagrid::make_point(trunk_size, trunk_size, total_height) * 1.01; }
};


struct BridgeConfig
{
  double beam_thickness;
  double beam_length;
  int instances;
};




struct BenchmarkResult
{
  BenchmarkResult() : time(0), volume_time_templated(0), facet_time_templated(0), line_time_templated(0), statistics_valid(false), statistics_valid_templated(false) {}

  double time;

  double volume_time_templated;
  double facet_time_templated;
  double line_time_templated;

  double time_templated() const { return line_time_templated + facet_time_templated + volume_time_templated; }

  shared_ptr<TemplatedMesh> templated_mesh;
  MeshType mesh;

  bool mr() const
  {
    RegionRangeType regions(mesh);
    return regions.size() > 1;
  }

  int cell_count() const { return viennagrid::cells(mesh).size(); }
  int size() const { return mr() ? mr_mesh_size(mesh) : mesh_size(mesh); }
  int size_CSR() const { return size() + mesh_matrix.memory_usage(); }

  int structure_instance_cell_count() const { return viennagrid::cells(templated_mesh->structure_instance).size(); }
  int size_templated() const { return templated_mesh->size(); }
  int size_templated_SVB() const { return templated_mesh->size_SVB(); }
  int size_templated_structure_instance() const { return mr() ? mr_mesh_size(templated_mesh->structure_instance) : mesh_size(templated_mesh->structure_instance); }

  int size_templated_CSR() const { return size_templated() + SI_matrix.memory_usage(); }
  int size_templated_SVB_CSR() const { return size_templated_SVB() + SI_matrix.memory_usage(); }
  int size_templated_structure_instance_CSR() const { return size_templated_structure_instance() + SI_matrix.memory_usage(); }
  int size_templated_CSR_templated() const { return size_templated() + templated_matrix.memory_usage(); }
  int size_templated_SVB_CSR_templated() const { return size_templated_SVB() + templated_matrix.memory_usage(); }

  double cell_count_ratio() const { return static_cast<double>(structure_instance_cell_count()) / static_cast<double>(cell_count()); }

  double benefit_time() const { return time/time_templated(); }
  double benefit_time_per_cell() const { return benefit_time() * cell_count_ratio(); }


  double structure_instance_benefit_size() const
  { return static_cast<double>(size_templated_structure_instance()) / static_cast<double>(size_templated()); }
  double structure_instance_benefit_size_SVB() const
  { return static_cast<double>(size_templated_structure_instance()) / static_cast<double>(size_templated_SVB()); }
  double structure_instance_benefit_size_CSR() const
  { return static_cast<double>(size_templated_structure_instance_CSR()) / static_cast<double>(size_templated_CSR()); }
  double structure_instance_benefit_size_SVB_CSR() const
  { return static_cast<double>(size_templated_structure_instance_CSR()) / static_cast<double>(size_templated_SVB_CSR()); }
  double structure_instance_benefit_size_CSR_templated() const
  { return static_cast<double>(size_templated_structure_instance_CSR()) / static_cast<double>(size_templated_CSR_templated()); }
  double structure_instance_benefit_size_SVB_CSR_templated() const
  { return static_cast<double>(size_templated_structure_instance_CSR()) / static_cast<double>(size_templated_SVB_CSR_templated()); }


  double remeshed_benefit_size() const
  { return static_cast<double>(size()) / static_cast<double>(size_templated()); }
  double remeshed_benefit_size_CSR() const
  { return static_cast<double>(size_CSR()) / static_cast<double>(size_templated_CSR()); }
  double remeshed_benefit_size_CSR_templated() const
  { return static_cast<double>(size_CSR()) / static_cast<double>(size_templated_CSR_templated()); }
  double remeshed_benefit_size_per_cell() const
  { return remeshed_benefit_size() * cell_count_ratio(); }
  double remeshed_benefit_size_per_cell_CSR() const
  { return remeshed_benefit_size_CSR() * cell_count_ratio(); }
  double remeshed_benefit_size_per_cell_CSR_templated() const
  { return remeshed_benefit_size_CSR_templated() * cell_count_ratio(); }

  double remeshed_benefit_size_SVB() const
  { return static_cast<double>(size()) / static_cast<double>(size_templated_SVB()); }
  double remeshed_benefit_size_SVB_CSR() const
  { return static_cast<double>(size_CSR()) / static_cast<double>(size_templated_SVB_CSR()); }
  double remeshed_benefit_size_SVB_CSR_templated() const
  { return static_cast<double>(size_CSR()) / static_cast<double>(size_templated_SVB_CSR_templated()); }
  double remeshed_benefit_size_SVB_per_cell() const
  { return remeshed_benefit_size_SVB() * cell_count_ratio(); }
  double remeshed_benefit_size_SVB_per_cell_CSR() const
  { return remeshed_benefit_size_SVB_CSR() * cell_count_ratio(); }
  double remeshed_benefit_size_SVB_per_cell_CSR_templated() const
  { return remeshed_benefit_size_SVB_CSR_templated() * cell_count_ratio(); }



  CSR SI_matrix;
  double time_SI_matrix;

  TCSR templated_matrix;
  double time_templated_matrix;

  CSR mesh_matrix;
  double time_mesh_matrix;

  void bench_matrix(int bench_count)
  {
    if (bench_count > 0)
    {
      viennautils::Timer timer;
      {
        std::vector<double> b = randomVector(SI_matrix.size());
        std::vector<double> result;
        timer.start();
        for (int i = 0; i != bench_count; ++i)
        {
          prod(SI_matrix, b, result);
        }
        time_SI_matrix = timer.get() / bench_count;
      }

      {
        std::vector<double> b = randomVector(templated_matrix.size);
        std::vector<double> result;
        timer.start();
        for (int i = 0; i != bench_count; ++i)
        {
          prod(templated_matrix, b, result);
        }
        time_templated_matrix = timer.get() / bench_count;
      }

      {
        std::vector<double> b = randomVector(mesh_matrix.size());
        std::vector<double> result;
        timer.start();
        for (int i = 0; i != bench_count; ++i)
        {
          prod(mesh_matrix, b, result);
        }
        time_mesh_matrix = timer.get() / bench_count;
      }
    }
  }




  bool statistics_valid;
  QuantityField angle_field;
  StatisticType angle_statistic;

  QuantityField ratio_field;
  StatisticType ratio_statistic;

  bool statistics_valid_templated;
  QuantityField angle_field_templated;
  StatisticType angle_statistic_templated;

  QuantityField ratio_field_templated;
  StatisticType ratio_statistic_templated;


  void finalize(int angle_hist_bin_count, int ratio_hist_bin_count);


  void save_conventional(std::string const & filename) const;
  void save_templated(std::string const & filename) const;
  void save(std::string const & filename) const;

  void print_conventional(BenchmarkConfig const & bench_cfg, std::ostream & stream = std::cout) const;
  void print_templated(BenchmarkConfig const & bench_cfg, std::ostream & stream = std::cout) const;
  void print_benefits(std::ostream & stream = std::cout) const;
  void print(BenchmarkConfig const & bench_cfg, std::ostream & stream = std::cout) const;
};


MeshType refine_lines(MeshType const & mesh, viennagrid_numeric line_size);

BenchmarkResult bench_reflection_2d_x(BenchmarkConfig const & bench_cfg, MeshType const & boundary_geometry);

BenchmarkResult bench_reflection_3d_xy(BenchmarkConfig const & bench_cfg, viennagrid_plc boundary_geometry);
BenchmarkResult bench_reflection_3d_xy(BenchmarkConfig const & bench_cfg, MeshType const & boundary_geometry);

BenchmarkResult bench_rotational_2d_circle(BenchmarkConfig const & bench_cfg, double radius);
BenchmarkResult bench_rotational_2d_circle_opt(BenchmarkConfig const & bench_cfg, double radius, double inner_radius);

BenchmarkResult bench_rotational_3d_tsv(BenchmarkConfig const & bench_cfg, TSVconfig const & cfg);
BenchmarkResult bench_rotational_3d_tsv_nonsym(BenchmarkConfig const & bench_cfg, TSVconfig const & cfg);

shared_ptr<TemplatedMesh> make_bridge_2d(BridgeConfig const & bridge_cfg, BenchmarkConfig const & bench_cfg,
                                         double * line_time, double * volume_time,
                                         viennamesh::point_container & hole_points);
BenchmarkResult bench_bridge_2d(BridgeConfig const & bridge_cfg, BenchmarkConfig const & bench_cfg);

BenchmarkResult bench_multi_tsv(BenchmarkConfig const & bench_cfg, TSVconfig const & cfg);


void bench_conventional_2d(BenchmarkConfig const & bench_cfg,
                           BenchmarkResult & bench_result,
                           MeshType const & structure_instance,
                           viennamesh::point_container const & hole_points = viennamesh::point_container());


#endif
