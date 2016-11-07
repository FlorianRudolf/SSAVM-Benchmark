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
  BenchmarkConfig bench_cfg(2);
  bench_cfg.from_args(argc, argv);

  viennamesh_log_set_info_level(-1);
  viennamesh_log_set_error_level(10);
  viennamesh_log_set_warning_level(-1);
  viennamesh_log_set_debug_level(-1);
  viennamesh_log_set_stack_level(-1);


  BridgeConfig bridge_cfg;
  bridge_cfg.beam_thickness = 10;
  bridge_cfg.beam_length = 100;
  bridge_cfg.instances = 7;



  {
    std::ofstream br_file( (bench_cfg.full_output_filename() + "_results").c_str() );
    bench_cfg.print();
    bench_cfg.print(br_file);

    BenchmarkResult br = bench_bridge_2d(bridge_cfg, bench_cfg);

    br.finalize(bench_cfg.angle_hist_bin_count, bench_cfg.ratio_hist_bin_count);
    br.bench_matrix(bench_cfg.matrix_bench_count);
    br.print(bench_cfg);
    br.print(bench_cfg, br_file);

    if (bench_cfg.output_active())
      br.save(bench_cfg.full_output_filename());
  }


  return 0;
}
