#include <string>
#include <sstream>

#include <tclap/CmdLine.h>
#include "pugixml.hpp"

#include "common.hpp"
#include "templated_mesh.hpp"
#include "benchmark.hpp"

#include "viennagrid/viennagrid.hpp"




int main(int argc, char **argv)
{
  BenchmarkConfig bench_cfg(2);
  bench_cfg.from_args(argc, argv);


  viennamesh_log_set_info_level(-1);
  viennamesh_log_set_error_level(10);
  viennamesh_log_set_warning_level(-1);
  viennamesh_log_set_debug_level(-1);
  viennamesh_log_set_stack_level(-1);

  std::ofstream br_file( (bench_cfg.full_output_filename() + "_results").c_str() );
  bench_cfg.print();
  bench_cfg.print(br_file);

  BenchmarkResult br = bench_rotational_2d_circle(bench_cfg, 1);
  br.finalize(bench_cfg.angle_hist_bin_count, bench_cfg.ratio_hist_bin_count);
  br.bench_matrix(bench_cfg.matrix_bench_count);
  br.print_conventional(bench_cfg);
  br.print_conventional(bench_cfg, br_file);
  br.print_templated(bench_cfg);
  br.print_templated(bench_cfg, br_file);

  if (bench_cfg.output_active())
    br.save(bench_cfg.full_output_filename());



  BenchmarkResult br_nonsym = bench_rotational_2d_circle_opt(bench_cfg, 1, bench_cfg.inner_radius());
  br_nonsym.finalize(bench_cfg.angle_hist_bin_count, bench_cfg.ratio_hist_bin_count);
  br_nonsym.bench_matrix(bench_cfg.matrix_bench_count);
  br_nonsym.print_templated(bench_cfg);
  br_nonsym.print_templated(bench_cfg, br_file);

  if (bench_cfg.output_active())
    br_nonsym.save_templated(bench_cfg.full_output_filename() + "_nonsym");

  br.print_benefits();
  br.print_benefits(br_file);
  br_nonsym.print_benefits();
  br_nonsym.print_benefits(br_file);

  return 0;
}
