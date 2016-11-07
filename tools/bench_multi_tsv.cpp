#include <string>
#include <sstream>

#include <tclap/CmdLine.h>
#include "pugixml.hpp"

#include "common.hpp"
#include "templated_mesh.hpp"
#include "benchmark.hpp"
#include "3d_geometry.hpp"






int main(int argc, char **argv)
{
  BenchmarkConfig bench_cfg(3);
  bench_cfg.from_args(argc, argv);

  viennamesh_log_set_info_level(-1);
  viennamesh_log_set_error_level(10);
  viennamesh_log_set_warning_level(-1);
  viennamesh_log_set_debug_level(-1);
  viennamesh_log_set_stack_level(-1);



  TSVconfig cfg;
  cfg.radius = 100;
  cfg.total_height = 270;
  cfg.hole_radius = 50;
  cfg.hole_height = 250;
  cfg.layer_thickness = 5;
  cfg.nonsym_inner_radius = 10;
  cfg.rotational_symmetry_order = bench_cfg.rotational_symmetry_order;
  cfg.distance = 2.2*cfg.radius;
  cfg.trunk_size = 1.5*cfg.distance;

  {
    std::ofstream br_file( (bench_cfg.full_output_filename() + "_results").c_str() );
    bench_cfg.print();
    bench_cfg.print(br_file);

    BenchmarkResult br = bench_multi_tsv(bench_cfg, cfg);

    br.finalize(bench_cfg.angle_hist_bin_count, bench_cfg.ratio_hist_bin_count);
    br.bench_matrix(bench_cfg.matrix_bench_count);
    br.print(bench_cfg);
    br.print(bench_cfg, br_file);

    if (bench_cfg.output_active())
      br.save(bench_cfg.full_output_filename());
  }

  return 0;
}
