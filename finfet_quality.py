#!/usr/bin/python

import sys
import os
import math
import Gnuplot as gp

from common import *



def check_quality( path ):

  #conventional_not_fulfilled = 0
  #templated_not_fulfilled = 0
  #both_not_fulfilled = 0
  #templated_worse = 0
  templated_better = 0
  templated_worst_better = 0
  total = 0

  print("\\hline")

  for file in os.listdir(path):

    if "_result" in file:

      filename = path + "/" + file

      config_cell_size = float(get_value(filename, "cell size"))
      config_min_angle = round(float(get_value(filename, "min angle")) * 180 / math.pi)
      config_radius_edge_ratio = float(get_value(filename, "radius edge ratio"))

      line = "$" + str(config_cell_size) + "$ & $"
      if config_min_angle >= 0:
        line += str(config_min_angle) + "\degree"
      else:
        line += str(config_radius_edge_ratio)
      line += "$ & "

      in_use = False

      if config_min_angle >= 0:
        total += 1

        worse_than_config = float(get_value(filename, "Under config value", 1))
        worse_than_config_templated = float(get_value(filename, "Under config value", 2))
        worst_quality = round(float(get_value(filename, "min", 2)) * 180 / math.pi, 3)
        worst_quality_templated = round(float(get_value(filename, "min", 4)) * 180 / math.pi, 3)

        relative_worst_value = round(worst_quality_templated/config_min_angle, 3)

        line += "$\\num{" + str(worse_than_config) + "}$ & $\\num{" + str(worse_than_config_templated) + "}$ & $" + str(worst_quality) + "\degree$ & $" + str(worst_quality_templated) + "\degree$ & $"

        if relative_worst_value > 1.3 or relative_worst_value < 0.7:
          line += "\\color{red}"
        line += str(relative_worst_value) + "$"

        if worse_than_config_templated > worse_than_config:
          in_use = True

        print()
        if worse_than_config_templated < worse_than_config:
          templated_better += 1
        if worst_quality < worst_quality_templated:
          templated_worst_better += 1


      if config_radius_edge_ratio != -1:
        worse_than_config = float(get_value(filename, "Over config value", 1))
        worse_than_config_templated = float(get_value(filename, "Over config value", 2))
        worst_quality = round(float(get_value(filename, "max", 2)), 3)
        worst_quality_templated = round(float(get_value(filename, "max", 4)), 3)

        relative_worst_value = round(worst_quality_templated/config_radius_edge_ratio, 3)

        line += "$\\num{" + str(worse_than_config) + "}$ & $\\num{" + str(worse_than_config_templated) + "}$ & $" + str(worst_quality) + "$ & $" + str(worst_quality_templated) + "$ & $"

        if relative_worst_value > 1.3 or relative_worst_value < 0.7:
          line += "\\color{red}"

        line += str(relative_worst_value) + "$"

        if worse_than_config_templated > worse_than_config:
          in_use = True

        #if worse_than_config_templated < worse_than_config:
          #templated_better += 1
        #if worst_quality < worst_quality_templated:
          #templated_worst_better += 1


      line += " \\\\"

      if in_use:
        print(line)
        print("\\hline")

  print("total = " + str(total) + "  templated better = " + str(templated_better) + "  templated worst better = " + str(templated_worst_better))


check_quality("finfet_3d")
