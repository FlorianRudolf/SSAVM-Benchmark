#!/usr/bin/python

import sys
import os
import Gnuplot as gp
import math

from common import *


pointstyle = "ps 2 lw 3"


def shift_pairs( pairs, factor ):
  for i in range(0, len(pairs)):
    pairs[i] = (factor*pairs[i][0], pairs[i][1])
  return pairs

def min_max_pairs(*args):
  l = list(args[0])
  for i in range(1, len(args)):
    l += args[i]

  min_ = min(l)
  max_ = max(l)

  return [ [min_[0], max_[0]], [min_[1], max_[1]] ]


def get_plot_data_cc( path, value_name, average = True, number = 1 ):
  pairs=[]
  for file in os.listdir(path):

    if "_result" in file:
      filename = path + "/" + file
      if not file_valid(filename):
        continue

      cc_value = float(get_value(filename, "cell count"))
      ccsi_value = float(get_value(filename, "Templated structure instance cell count"))

      second_value = float(get_value(filename, value_name, number))

      if average:
        pairs.append( ((cc_value+ccsi_value)/2, second_value) )
      else:
        pairs.append( (ccsi_value, second_value) )

  #pairs = sorted(pairs, key=lambda tup: tup[0])
  return pairs

def get_plot_data_rso( path, value_name, number = 1, cc_threshold = 0 ):
  pairs=[]
  for file in os.listdir(path):

    if "_result" in file:
      filename = path + "/" + file
      if not file_valid(filename):
        continue

      cc = (float(get_value(filename, "cell count")) + float(get_value(filename, "Templated structure instance cell count")) + float(get_value(filename, "Templated structure instance cell count", 2))) / 3

      if cc >= cc_threshold:
        rso_value = float(get_value(filename, "rotational symmetry order"))
        second_value = float(get_value(filename, value_name, number))
        pairs.append( (rso_value, second_value) )

  return pairs



def get_histogram( filename, number = 1 ):
  f = open(filename, 'r')
  i = 1
  go = False

  pairs = []

  for line in f:
    if go:
      if line.find("  [") != 0:
        return sorted(pairs, key=lambda tup: tup[0])

      tokens = line.split()
      value = tokens[-1]
      quality = tokens[0]

      left = quality[1:quality.find(",")]
      right = quality[quality.find(",")+1:quality.find("]")]

      if left == "-inf":
        continue

      pairs.append( (float(left) * 180/math.pi, float(value)) )

    if line.find("median") != 0: continue
    if i == number:
      go = True
    else:
      i = i+1


def plot_histogram( filename, number, output_filename, good_bin, histogram_title = "", yrange = [0.001, 0.3], size_x = 1, size_y = 0.6 ):
  g = gp.Gnuplot()
  g("set terminal postscript eps enhanced color font 'Helvetica,36'")
  g("set output '" + output_filename + "'")
  g("set size " + str(size_x) + "," + str(size_y))

  g("set style data histogram")
  g("set style histogram")
  g("set style fill solid 1.0 border -1")
  if histogram_title != "":
    g("set title \"" + histogram_title + "\"")

  data = get_histogram(filename, number)
  b_good = gp.Data(data[0:good_bin], with_="boxes lc rgb 'red' lt 1", title = "")
  b_bad = gp.Data(data[good_bin:], with_="boxes lc rgb 'green' lt 1", title = "")

  x = (data[good_bin-1][0] + data[good_bin][0])/2

  g("set arrow from " + str(x) + "," + str(yrange[0]) + " to " + str(x) + ",0.3 nohead")
  g("set logscale y")
  g("set xlabel \"Smallest angle\"")
  #g("set xrange [" + str(data[0][0]) + ":" + str(data[-1][0]) + "]")
  g("set xrange [0:60]")
  #g("set ytics (\"0\" 0.01, \"1\" 1, \"10\" 10)")

  g("set yrange [" + str(yrange[0]) + ":" + str(yrange[1]) + "]")

  g.plot(b_good, b_bad)

  return g

def plot_hists( path, output_path ):

  f = open(output_path + "/latex_figures.tex", "w")
  latex_path = output_path[output_path.find("document/")+9:]
  #print(latex_path)

  #f.write("\\begin{figure}[t]")
  #f.write("")

  number = 0
  for file in os.listdir(path):

    if "_result" in file:
      filename = path + "/" + file
      if not file_valid(filename):
        continue

      cell_size = float(get_value(filename, "cell size"))
      angle = round(float(get_value(filename, "min angle")) / math.pi * 180)
      rso = int(get_value(filename, "rotational symmetry order"))

      under_str = get_value(filename, "Under config value", 2)
      if not under_str:
        continue

      under = float(under_str)

      #if angle < 30  or rso <= 64:
        #continue

      if under <= 0.015:
        continue

      number += 1
      #if number == 10:
        #return

      good_bin = int(angle / 90 * 18)
      plot_histogram( filename, 1, output_path + "/" + str(number) + "_conventional.eps", good_bin, "Conventional", [0.001, 0.3], 1, 0.8 )
      plot_histogram( filename, 3, output_path + "/" + str(number) + "_templated.eps", good_bin, "Structure instance", [0.001, 0.3], 1, 0.8 )
      plot_histogram( filename, 5, output_path + "/" + str(number) + "_templated_sao.eps", good_bin, "SAO", [0.001, 0.3], 1, 0.8 )

      f.write("  \\begin{subfigure}[b]{0.95\\textwidth}\n")
      f.write("    \\begin{center}\n")
      f.write("      \\centering\n")
      f.write("      \\includegraphics[width=0.32\\textwidth]{" + latex_path + "/" + str(number) + "_conventional}\n")
      f.write("      \\includegraphics[width=0.32\\textwidth]{" + latex_path + "/" + str(number) + "_templated}\n")
      f.write("      \\includegraphics[width=0.32\\textwidth]{" + latex_path + "/" + str(number) + "_templated_sao}\n")
      f.write("    \\end{center}\n")
      f.write("    \\vspace{-0.3cm}\n")
      f.write("    \\caption{$" + str(rso) + "$-polygon, smallest angle = $" + str(angle) + "\\degree$, cell size = $\\num{" + str(cell_size) + "}$}\n")
      f.write("  \\end{subfigure}\n")

  f.close()




def plot_1( path, expected_benefit, output_filename, name = "benefit time", number = 1, keypos = "right bottom" ):
  g = gp.Gnuplot()
  g("set terminal postscript eps enhanced color font 'Helvetica,24'")
  g("set output '" + output_filename + "'")

  data = get_plot_data_cc(path, name, number)

  min_max = min_max_pairs(data)
  min_max[0][0] *= 0.9
  min_max[0][1] *= 1.1
  min_max[1][0] -= 2*expected_benefit/10
  min_max[1][1] += 2*expected_benefit/10

  b = gp.Data(data, with_="points lc rgb 'red' " + pointstyle + " pt 2", title = "")

  g("set xrange [" + str(min_max[0][0]) + ":" + str(min_max[0][1]) + "]")
  #g("set yrange [" + str(min_y) + ":" + str(max_y) + "]")
  g("set logscale x")
  g("set format x \"10^{%L}\"")
  g("set key " + keypos)

  g("set xlabel \"Cell count\"")
  g("set ylabel \"Speedup\"")

  if expected_benefit > 0:
    g.plot(str(expected_benefit) + " with lines lc rgb 'green' notitle", b)
  else:
    g.plot(b)

  return g


def plot_2( path, expected_benefit, output_filename, name_1 = "benefit remeshed size", number_1 = 1, title_1 = "w/o SVB", name_2 = "benefit remeshed size SVB", number_2 = 1, title_2 = "w/ SVB", average = True, keypos = "right bottom", y_range = [-1,-1]):
  g = gp.Gnuplot()
  g("set terminal postscript eps enhanced color font 'Helvetica,24'")
  g("set output '" + output_filename + "'")

  data = get_plot_data_cc(path, name_1, average, number_1)
  data_SVB = get_plot_data_cc(path, name_2, average, number_2)

  min_max = min_max_pairs(data, data_SVB)

  min_max[0][0] *= 0.9
  min_max[0][1] *= 1.1
  min_max[1][0] -= 2*expected_benefit/10
  min_max[1][1] += 2*expected_benefit/10

  b = gp.Data(data, with_="points lc rgb 'red' " + pointstyle + " pt 8", title = title_1)
  b_SVB = gp.Data(data_SVB, with_="points lc rgb 'blue' " + pointstyle + " pt 6", title = title_2)

  g("set xrange [" + str(min_max[0][0]) + ":" + str(min_max[0][1]) + "]")
  if y_range != [-1,-1]:
    g("set yrange [" + str(y_range[0]) + ":" + str(y_range[1]) + "]")
  g("set logscale x")
  g("set format x \"10^{%L}\"")
  g("set key " + keypos)

  g("set xlabel \"Cell count\"")
  g("set ylabel \"Memory savings\"")

  if expected_benefit > 0:
    g.plot(str(expected_benefit) + " with lines lc rgb 'green' notitle", b, b_SVB)
  else:
    g.plot(b, b_SVB)

  return g




def plot_time_benefits_rotational_rso( path, output_filename, cc_threshold, color, color_mao ):
  g = gp.Gnuplot()
  g("set terminal postscript eps enhanced color font 'Helvetica,24'")
  g("set output '" + output_filename + "'")

  data = shift_pairs(get_plot_data_rso(path, "benefit time", 1, cc_threshold), 0.96)
  data_nonsym = shift_pairs(get_plot_data_rso(path, "benefit time", 2, cc_threshold), 1.04)

  min_max = min_max_pairs(data, data_nonsym)
  min_max[0][0] *= 0.9
  min_max[0][1] *= 1.1

  b = gp.Data(data, with_="points lc rgb '" + color + "' " + pointstyle + " pt 2", title = "templated slice")
  b_nonsym = gp.Data(data_nonsym, with_="points lc rgb '" + color_mao + "' " + pointstyle + " pt 2", title = "small angle optimization")

  g("set xrange [" + str(min_max[0][0]) + ":" + str(min_max[0][1]) + "]")
  #g("set yrange [" + str(min_y) + ":" + str(max_y) + "]")
  g("set logscale x")
  g("set logscale y")
  g("set format x \"10^{%L}\"")
  g("set key top left")

  g("set xlabel \"Rotational Symmetry Order\"")
  g("set ylabel \"Speedup\"")

  #g("set style data points")
  g.plot("x with lines lc rgb 'black' lt 1 notitle", "1 with lines lc rgb 'black' lt 1 notitle",
         b, b_nonsym)

  return g


def plot_size_benefits_rotational_rso( path, output_filename, cc_threshold, color, color_SVB, color_mao, color_SVB_mao ):
  g = gp.Gnuplot()
  g("set terminal postscript eps enhanced color font 'Helvetica,24'")
  g("set output '" + output_filename + "'")

  data = shift_pairs(get_plot_data_rso(path, "benefit remeshed size", 1, cc_threshold), 0.88)
  data_SVB = shift_pairs(get_plot_data_rso(path, "benefit remeshed size SVB", 1, cc_threshold), 0.96)

  data_nonsym = shift_pairs(get_plot_data_rso(path, "benefit remeshed size", 2, cc_threshold), 1.04)
  data_SVB_nonsym = shift_pairs(get_plot_data_rso(path, "benefit remeshed size SVB", 2, cc_threshold), 1.12)

  min_max = min_max_pairs(data, data_SVB, data_nonsym, data_SVB_nonsym)
  min_max[0][0] *= 0.9
  min_max[0][1] *= 1.1

  b = gp.Data(data, with_="points lc rgb '" + color + "' " + pointstyle + " pt 8", title = "TS w/o SVB")
  b_SVB= gp.Data(data_SVB, with_="points lc rgb '" + color_SVB + "' " + pointstyle + " pt 6", title = "TS w/ SVB")
  b_nonsym = gp.Data(data_nonsym, with_="points lc rgb '" + color_mao + "' " + pointstyle + " pt 8", title = "SAO w/o SVB")
  b_SVB_nonsym = gp.Data(data_SVB_nonsym, with_="points lc rgb '" + color_SVB_mao + "' " + pointstyle + " pt 6", title = "SAO w/ SVB")

  g("set xrange [" + str(min_max[0][0]) + ":" + str(min_max[0][1]) + "]")
  #g("set yrange [" + str(min_y) + ":" + str(max_y) + "]")
  g("set logscale x")
  g("set logscale y")
  g("set format x \"10^{%L}\"")
  g("set key top left")

  g("set xlabel \"Rotational Symmetry Order\"")
  g("set ylabel \"Memory savings\"")

  #g("set style data points")
  g.plot("x with lines lc rgb 'black' lt 1 notitle", "1 with lines lc rgb 'black' lt 1 notitle",
         b, b_SVB, b_nonsym, b_SVB_nonsym)

  return g





def plot_time_benefits_rotational_cc( path, output_filename, number, color ):
  g = gp.Gnuplot()
  g("set terminal postscript eps enhanced color font 'Helvetica,24'")
  g("set output '" + output_filename + "'")

  data = get_plot_data_cc(path, "benefit time", True, number)
  #data_nonsym = get_plot_data_cc(path, "benefit time", True, 2)


  min_max = min_max_pairs(data)
  min_max[0][0] *= 0.9
  min_max[0][1] *= 1.1

  b = gp.Data(data, with_="points lc rgb '" + color + "' " + pointstyle + " pt 2", title = "")
  #b_nonsym = gp.Data(data_nonsym, with_="points lc rgb 'blue' ps 2 pt 2", title = "")

  g("set xrange [" + str(min_max[0][0]) + ":" + str(min_max[0][1]) + "]")
  g("set yrange [0.1:1000]")
  g("set logscale x")
  g("set logscale y")
  g("set format x \"10^{%L}\"")
  g("set key right bottom")

  g("set xlabel \"Cell count\"")
  g("set ylabel \"Speedup\"")

  #g("set style data points")
  g.plot("8 with lines lc rgb 'black' lt 2 notitle", "16 with lines lc rgb 'black' lt 2 notitle", "32 with lines lc rgb 'black' lt 2 notitle",
         "64 with lines lc rgb 'black' lt 2 notitle", "128 with lines lc rgb 'black' lt 2 notitle", "1 with lines lc rgb 'black' lt 1 notitle",
         b)

  return g


def plot_size_benefits_rotational_cc( path, output_filename, number,
                                     color1, color2,
                                     value_name1 = "benefit remeshed size", value_name2 = "benefit remeshed size SVB",
                                     title1 = "w/o SVB", title2 = "w/ SVB"):
  g = gp.Gnuplot()
  g("set terminal postscript eps enhanced color font 'Helvetica,24'")
  g("set output '" + output_filename + "'")

  data = get_plot_data_cc(path, value_name1, True, number)
  data_SVB = get_plot_data_cc(path, value_name2, True, number)

  #data_nonsym = get_plot_data_cc(path, "benefit remeshed size", True, 2)
  #data_SVB_nonsym = get_plot_data_cc(path, "benefit remeshed size SVB", True, 2)

  min_max = min_max_pairs(data, data_SVB)
  min_max[0][0] *= 0.9
  min_max[0][1] *= 1.1

  b = gp.Data(data, with_="points lc rgb '" + color1 + "' " + pointstyle + " pt 8", title = title1)
  b_SVB= gp.Data(data_SVB, with_="points lc rgb '" + color2 + "' " + pointstyle + " pt 6", title = title2)
  #b_nonsym = gp.Data(data_nonsym, with_="points lc rgb 'blue' ps 2 pt 2", title = "")
  #b_SVB_nonsym = gp.Data(data_SVB_nonsym, with_="points lc rgb 'blue' ps 2 pt 6", title = "")

  g("set xrange [" + str(min_max[0][0]) + ":" + str(min_max[0][1]) + "]")
  g("set yrange [0.1:1000]")
  g("set logscale x")
  g("set logscale y")
  g("set format x \"10^{%L}\"")
  g("set key right bottom")

  g("set xlabel \"Cell count\"")
  g("set ylabel \"Memory savings\"")

  #g("set style data points")
  g.plot("8 with lines lc rgb 'black' lt 2 notitle", "16 with lines lc rgb 'black' lt 2 notitle", "32 with lines lc rgb 'black' lt 2 notitle",
         "64 with lines lc rgb 'black' lt 2 notitle", "128 with lines lc rgb 'black' lt 2 notitle", "1 with lines lc rgb 'black' lt 1 notitle",
         b, b_SVB)

  return g









plot_1("aircraft", 2, "../document/figures/benchmark/aircraft_time.eps")
plot_2("aircraft", 2, "../document/figures/benchmark/aircraft_size.eps", "benefit remeshed size", 1, "w/o SVB", "benefit remeshed size SVB", 1, "w/ SVB", True, "right bottom", [0.8,2.1])
plot_2("aircraft", 2, "../document/figures/benchmark/aircraft_size_SI.eps", "benefit SI size", 1, "w\o SVB", "benefit SI size SVB", 1, "w\ SVB", False, "right bottom", [1.4,2.1])

plot_1("mosfet_2d", 2, "../document/figures/benchmark/mosfet_2d_time.eps")
plot_2("mosfet_2d", 2, "../document/figures/benchmark/mosfet_2d_size.eps")
plot_2("mosfet_2d", 2, "../document/figures/benchmark/mosfet_2d_size_SI.eps", "benefit SI size", 1, "w\o SVB", "benefit SI size SVB", 1, "w\ SVB", False)

plot_1("finfet_3d", 4, "../document/figures/benchmark/finfet_3d_time.eps")
plot_2("finfet_3d", 4, "../document/figures/benchmark/finfet_3d_size.eps")
plot_2("finfet_3d", 4, "../document/figures/benchmark/finfet_3d_size_SI.eps", "benefit SI size", 1, "w\o SVB",  "benefit SI size SVB", 1, "w\ SVB", False)

plot_1("bridge", -1, "../document/figures/benchmark/bridge_time.eps")
plot_2("bridge", -1, "../document/figures/benchmark/bridge_size.eps")
#plot_2("bridge", 4, "../document/figures/benchmark/bridge_size_SI.eps", "benefit SI size", 1, "w\o SVB",  "benefit SI size SVB", 1, "w\ SVB", False)

plot_1("multi_tsv", -1, "../document/figures/benchmark/multi_tsv_time.eps")
plot_2("multi_tsv", -1, "../document/figures/benchmark/multi_tsv_size.eps")
#plot_2("multi_tsv", 4, "../document/figures/benchmark/multi_tsv_size_SI.eps", "benefit SI size", 1, "w\o SVB",  "benefit SI size SVB", 1, "w\ SVB", False)


plot_time_benefits_rotational_rso("circle", "../document/figures/benchmark/circle_time_rso.eps", 100000, "red", "purple")
plot_size_benefits_rotational_rso("circle", "../document/figures/benchmark/circle_size_rso.eps", 5000, "red", "blue", "purple", "green")
plot_time_benefits_rotational_cc("circle", "../document/figures/benchmark/circle_time_cc.eps", 1, "red")
plot_size_benefits_rotational_cc("circle", "../document/figures/benchmark/circle_size_cc.eps", 1, "red", "blue")
plot_time_benefits_rotational_cc("circle", "../document/figures/benchmark/circle_time_cc_sao.eps", 2, "purple")
plot_size_benefits_rotational_cc("circle", "../document/figures/benchmark/circle_size_cc_sao.eps", 2, "purple", "green")



plot_time_benefits_rotational_rso("tsv", "../document/figures/benchmark/tsv_time_rso.eps", 0, "red", "purple")
plot_size_benefits_rotational_rso("tsv", "../document/figures/benchmark/tsv_size_rso.eps", 0, "red", "blue", "purple", "green")
plot_time_benefits_rotational_cc("tsv", "../document/figures/benchmark/tsv_time_cc.eps", 1, "red")
plot_size_benefits_rotational_cc("tsv", "../document/figures/benchmark/tsv_size_cc.eps", 1, "red", "blue")
plot_time_benefits_rotational_cc("tsv", "../document/figures/benchmark/tsv_time_cc_sao.eps", 2, "purple")
plot_size_benefits_rotational_cc("tsv", "../document/figures/benchmark/tsv_size_cc_sao.eps", 2, "purple", "green")

plot_size_benefits_rotational_cc("circle/quality", "../document/figures/benchmark/circle_size_cc_matrix.eps", 1, "red", "blue", "benefit remeshed size", "CSR benefit remeshed size", "w/o system matrix", "w/ system matrix")
plot_size_benefits_rotational_cc("circle/quality", "../document/figures/benchmark/circle_size_cc_matrix_SVB.eps", 1, "red", "blue", "benefit remeshed size SVB", "CSR benefit remeshed size SVB", "w/o system matrix", "w/ system matrix")

plot_size_benefits_rotational_cc("tsv", "../document/figures/benchmark/tsv_size_cc_matrix.eps", 1, "red", "blue", "benefit remeshed size", "CSR benefit remeshed size", "w/o system matrix", "w/ system matrix")
plot_size_benefits_rotational_cc("tsv", "../document/figures/benchmark/tsv_size_cc_matrix_SVB.eps", 1, "red", "blue", "benefit remeshed size SVB", "CSR benefit remeshed size SVB", "w/o system matrix", "w/ system matrix")


plot_size_benefits_rotational_cc("circle/quality", "../document/figures/benchmark/circle_size_cc_matrix_templated.eps", 1, "red", "blue", "benefit remeshed size", "templated CSR benefit remeshed size", "w/o system matrix", "w/ templated system matrix")
plot_size_benefits_rotational_cc("circle/quality", "../document/figures/benchmark/circle_size_cc_matrix_SVB_templated.eps", 1, "red", "blue", "benefit remeshed size SVB", "templated CSR benefit remeshed size SVB", "w/o system matrix", "w/ templated system matrix")

plot_size_benefits_rotational_cc("tsv", "../document/figures/benchmark/tsv_size_cc_matrix_templated.eps", 1, "red", "blue", "benefit remeshed size", "templated CSR benefit remeshed size", "w/o system matrix", "w/ templated system matrix")
plot_size_benefits_rotational_cc("tsv", "../document/figures/benchmark/tsv_size_cc_matrix_SVB_templated.eps", 1, "red", "blue", "benefit remeshed size SVB", "templated CSR benefit remeshed size SVB", "w/o system matrix", "w/ templated system matrix")


plot_2("circle/quality", -1, "../document/figures/benchmark/circle_matrix_benchmarks.eps", "matrix-vector product templated CSR vs remeshed CSR", 1, "TS", "matrix-vector product templated CSR vs remeshed CSR", 2, "SAO", True, "left top")
plot_2("tsv/quality", -1, "../document/figures/benchmark/tsv_matrix_benchmarks.eps", "matrix-vector product templated CSR vs remeshed CSR", 1, "TS", "matrix-vector product templated CSR vs remeshed CSR", 2, "SAO", True, "left top")


plot_hists("circle/quality", "../document/figures/benchmark/circle_histograms")

plot_histogram( "aircraft/vis/aircraft_c1024_a0.523599_results", 1, "../document/figures/benchmark/aircraft_hist.eps", 6, "", [0.0008, 0.3] )
plot_histogram( "aircraft/vis/aircraft_c1024_a0.523599_results", 3, "../document/figures/benchmark/aircraft_hist_symmetric.eps", 6, "", [0.0008, 0.3] )
