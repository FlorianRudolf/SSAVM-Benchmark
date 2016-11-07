#!/usr/bin/python

import sys
import os
import Gnuplot as gp
import math

from common import *


def get_best( path, value_name, rso = -1, number = 1 ):
  best = -1

  for file in os.listdir(path):
    if "_result" in file:
      filename = path + "/" + file
      if not file_valid(filename):
        continue

      if rso != -1 and int(get_value(filename, "rotational symmetry order")) != rso:
        continue

      value = float(get_value(filename, value_name, number))

      if value > best:
        best = value

  #if rso > 0:
    #print(path + " " + str(rso) + " " + value_name + " " + str(best))
  #else:
    #print(path + " " + value_name + " " + str(best))
  return best

def print_best( path, rso = -1, number = 1 ):
  best_mem = get_best(path, "benefit remeshed size", rso)
  best_mem_SVB = get_best(path, "benefit remeshed size SVB", rso)
  best_time = get_best(path, "benefit time", rso)

  os = path

  if rso > 0:
    os += " " + str(rso)

  os += "    $" + str(best_mem) + "$ & $" + str(best_mem_SVB) + "$ & $" + str(best_time) + "$ \\\\"
  print(os)

#"benefit time"
#"benefit remeshed size"
#"benefit remeshed size SVB"


print_best("aircraft")
print_best("mosfet_2d")
print_best("finfet_3d")

for rso in [8, 16, 32, 64, 128]:
  print_best("circle", rso)

for rso in [8, 16, 32, 64, 128]:
  print_best("tsv", rso)

print_best("bridge")
print_best("multi_tsv")


#plot_1("aircraft", 2, "../document/figures/benchmark/aircraft_time.eps")
#plot_2("aircraft", 2, "../document/figures/benchmark/aircraft_size.eps", "benefit remeshed size", 1, "w/o SVB", "benefit remeshed size SVB", 1, "w/ SVB", True, "right bottom", [0.8,2.1])
#plot_2("aircraft", 2, "../document/figures/benchmark/aircraft_size_SI.eps", "benefit SI size", 1, "w\o SVB", "benefit SI size SVB", 1, "w\ SVB", False, "right bottom", [1.4,2.1])

#plot_1("mosfet_2d", 2, "../document/figures/benchmark/mosfet_2d_time.eps")
#plot_2("mosfet_2d", 2, "../document/figures/benchmark/mosfet_2d_size.eps")
#plot_2("mosfet_2d", 2, "../document/figures/benchmark/mosfet_2d_size_SI.eps", "benefit SI size", 1, "w\o SVB", "benefit SI size SVB", 1, "w\ SVB", False)

#plot_1("finfet_3d", 4, "../document/figures/benchmark/finfet_3d_time.eps")
#plot_2("finfet_3d", 4, "../document/figures/benchmark/finfet_3d_size.eps")
#plot_2("finfet_3d", 4, "../document/figures/benchmark/finfet_3d_size_SI.eps", "benefit SI size", 1, "w\o SVB",  "benefit SI size SVB", 1, "w\ SVB", False)

#plot_1("bridge", -1, "../document/figures/benchmark/bridge_time.eps")
#plot_2("bridge", -1, "../document/figures/benchmark/bridge_size.eps")
##plot_2("bridge", 4, "../document/figures/benchmark/bridge_size_SI.eps", "benefit SI size", 1, "w\o SVB",  "benefit SI size SVB", 1, "w\ SVB", False)

#plot_1("multi_tsv", -1, "../document/figures/benchmark/multi_tsv_time.eps")
#plot_2("multi_tsv", -1, "../document/figures/benchmark/multi_tsv_size.eps")
##plot_2("multi_tsv", 4, "../document/figures/benchmark/multi_tsv_size_SI.eps", "benefit SI size", 1, "w\o SVB",  "benefit SI size SVB", 1, "w\ SVB", False)


#plot_time_benefits_rotational_rso("circle", "../document/figures/benchmark/circle_time_rso.eps", 100000, "red", "purple")
#plot_size_benefits_rotational_rso("circle", "../document/figures/benchmark/circle_size_rso.eps", 5000, "red", "blue", "purple", "green")
#plot_time_benefits_rotational_cc("circle", "../document/figures/benchmark/circle_time_cc.eps", 1, "red")
#plot_size_benefits_rotational_cc("circle", "../document/figures/benchmark/circle_size_cc.eps", 1, "red", "blue")
#plot_time_benefits_rotational_cc("circle", "../document/figures/benchmark/circle_time_cc_sao.eps", 2, "purple")
#plot_size_benefits_rotational_cc("circle", "../document/figures/benchmark/circle_size_cc_sao.eps", 2, "purple", "green")



#plot_time_benefits_rotational_rso("tsv", "../document/figures/benchmark/tsv_time_rso.eps", 0, "red", "purple")
#plot_size_benefits_rotational_rso("tsv", "../document/figures/benchmark/tsv_size_rso.eps", 0, "red", "blue", "purple", "green")
#plot_time_benefits_rotational_cc("tsv", "../document/figures/benchmark/tsv_time_cc.eps", 1, "red")
#plot_size_benefits_rotational_cc("tsv", "../document/figures/benchmark/tsv_size_cc.eps", 1, "red", "blue")
#plot_time_benefits_rotational_cc("tsv", "../document/figures/benchmark/tsv_time_cc_sao.eps", 2, "purple")
#plot_size_benefits_rotational_cc("tsv", "../document/figures/benchmark/tsv_size_cc_sao.eps", 2, "purple", "green")

#plot_size_benefits_rotational_cc("circle/quality", "../document/figures/benchmark/circle_size_cc_matrix.eps", 1, "red", "blue", "benefit remeshed size", "CSR benefit remeshed size", "w/o system matrix", "w/ system matrix")
#plot_size_benefits_rotational_cc("circle/quality", "../document/figures/benchmark/circle_size_cc_matrix_SVB.eps", 1, "red", "blue", "benefit remeshed size SVB", "CSR benefit remeshed size SVB", "w/o system matrix", "w/ system matrix")

#plot_size_benefits_rotational_cc("tsv", "../document/figures/benchmark/tsv_size_cc_matrix.eps", 1, "red", "blue", "benefit remeshed size", "CSR benefit remeshed size", "w/o system matrix", "w/ system matrix")
#plot_size_benefits_rotational_cc("tsv", "../document/figures/benchmark/tsv_size_cc_matrix_SVB.eps", 1, "red", "blue", "benefit remeshed size SVB", "CSR benefit remeshed size SVB", "w/o system matrix", "w/ system matrix")


#plot_size_benefits_rotational_cc("circle/quality", "../document/figures/benchmark/circle_size_cc_matrix_templated.eps", 1, "red", "blue", "benefit remeshed size", "templated CSR benefit remeshed size", "w/o system matrix", "w/ templated system matrix")
#plot_size_benefits_rotational_cc("circle/quality", "../document/figures/benchmark/circle_size_cc_matrix_SVB_templated.eps", 1, "red", "blue", "benefit remeshed size SVB", "templated CSR benefit remeshed size SVB", "w/o system matrix", "w/ templated system matrix")

#plot_size_benefits_rotational_cc("tsv", "../document/figures/benchmark/tsv_size_cc_matrix_templated.eps", 1, "red", "blue", "benefit remeshed size", "templated CSR benefit remeshed size", "w/o system matrix", "w/ templated system matrix")
#plot_size_benefits_rotational_cc("tsv", "../document/figures/benchmark/tsv_size_cc_matrix_SVB_templated.eps", 1, "red", "blue", "benefit remeshed size SVB", "templated CSR benefit remeshed size SVB", "w/o system matrix", "w/ templated system matrix")


#plot_2("circle/quality", -1, "../document/figures/benchmark/circle_matrix_benchmarks.eps", "matrix-vector product templated CSR vs remeshed CSR", 1, "TS", "matrix-vector product templated CSR vs remeshed CSR", 2, "SAO", True, "left top")
#plot_2("tsv/quality", -1, "../document/figures/benchmark/tsv_matrix_benchmarks.eps", "matrix-vector product templated CSR vs remeshed CSR", 1, "TS", "matrix-vector product templated CSR vs remeshed CSR", 2, "SAO", True, "left top")


#plot_hists("circle/quality", "../document/figures/benchmark/circle_histograms")
