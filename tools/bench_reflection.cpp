#include <string>
#include <sstream>

#include "common.hpp"
#include "benchmark.hpp"

#include "viennameshpp/core.hpp"
#include "viennameshpp/algorithm_pipeline.hpp"


int main(int argc, char **argv)
{
  BenchmarkConfig bench_cfg(-1);
  bench_cfg.from_args(argc, argv);

  viennamesh_log_set_info_level(-1);
  viennamesh_log_set_error_level(-1);
  viennamesh_log_set_warning_level(-1);
  viennamesh_log_set_debug_level(-1);
  viennamesh_log_set_stack_level(-1);

  viennagrid::mesh template_geometry;
  viennagrid_plc template_plc = 0;

  if (bench_cfg.template_filename.find(".poly") != std::string::npos)
  {
    viennagrid_plc_create(&template_plc);
    viennagrid_plc_read_tetgen_poly(template_plc, bench_cfg.template_filename.c_str());
    bench_cfg.dim = 3;
  }
  else
  {
    viennamesh::algorithm_handle mesh_reader = context().make_algorithm("mesh_reader");
    mesh_reader.set_input("filename", bench_cfg.template_filename);
    mesh_reader.run();
    template_geometry = mesh_reader.get_output<viennagrid_mesh>("mesh")();
    bench_cfg.dim = 2;
  }


  std::ofstream br_file( (bench_cfg.full_output_filename() + "_results").c_str() );
  bench_cfg.print();
  bench_cfg.print(br_file);


  BenchmarkResult br;

  if (template_plc)
  {
    br = bench_reflection_3d_xy(bench_cfg, template_plc);
  }
  else
  {
    br = bench_reflection_2d_x(bench_cfg, template_geometry);
  }

  br.finalize(bench_cfg.angle_hist_bin_count, bench_cfg.ratio_hist_bin_count);
  br.bench_matrix(bench_cfg.matrix_bench_count);

  br.print(bench_cfg);
  br.print(bench_cfg, br_file);

  br_file.close();

  if (bench_cfg.output_active())
  {
    br.save(bench_cfg.full_output_filename());
  }

  if (template_plc)
    viennagrid_plc_release(template_plc);


  return 0;
}
