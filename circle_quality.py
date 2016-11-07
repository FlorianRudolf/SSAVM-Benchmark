#!/usr/bin/python

import sys
import os
import math
import Gnuplot as gp

from common import *



def check_quality( path ):


  for rso in [8,16,32,64,128]:

    total = 0
    conventional_not_fulfilled = 0
    templated_not_fulfilled = 0
    templated_nonsym_not_fulfilled = 0
    templated_nonsym_worse_than_templated = 0


    print("RSO = " + str(rso))

    for file in os.listdir(path):

      if "_result" in file:
        filename = path + "/" + file
        if not file_valid(filename):
          continue

        rotational_symmetry_order = float(get_value(filename, "rotational symmetry order"))
        if rotational_symmetry_order != rso:
          continue

        total += 1

        config_cell_size = float(get_value(filename, "cell size"))
        config_min_angle = round(float(get_value(filename, "min angle")) * 180 / math.pi)

        worse_than_config = float(get_value(filename, "Under config value", 1))
        worse_than_config_templated = float(get_value(filename, "Under config value", 2))
        worse_than_config_nonsym = float(get_value(filename, "Under config value", 3))
        worst_quality = round(float(get_value(filename, "min", 2)) * 180 / math.pi, 3)
        worst_quality_templated = round(float(get_value(filename, "min", 4)) * 180 / math.pi, 3)
        worst_quality_templated_nonsym = round(float(get_value(filename, "min", 6)) * 180 / math.pi, 3)

        if worse_than_config > 0:
          conventional_not_fulfilled += 1
        if worse_than_config_templated > 0:
          templated_not_fulfilled += 1
        if worse_than_config_nonsym > 0:
          templated_nonsym_not_fulfilled += 1
        if worse_than_config_templated < worse_than_config_nonsym:
          templated_nonsym_worse_than_templated += 1



        if worse_than_config_nonsym > 0:
          print(filename)
          print("  config = " + str(config_min_angle))
          if worse_than_config_templated < worse_than_config_nonsym:
            print("     !!!!!!!")
          print("  worse than config: " + str(worse_than_config) + "  templated: " + str(worse_than_config_templated) + "  nonsym: " + str(worse_than_config_nonsym))
          print("  worset quality: " + str(worst_quality) + "  templated: " + str(worst_quality_templated) + "  nonsym: " + str(worst_quality_templated_nonsym) + "\n")

    print("total = " + str(total))
    print("conventional_not_fulfilled = " + str(conventional_not_fulfilled))
    print("templated_not_fulfilled = " + str(templated_not_fulfilled))
    print("templated_nonsym_not_fulfilled = " + str(templated_nonsym_not_fulfilled))
    print("templated_nonsym_worse_than_templated = " + str(templated_nonsym_worse_than_templated))
    print("\n\n")


        #relative_worst_value = round(worst_quality_templated/config_min_angle, 3)

        #line += "$\\num{" + str(worse_than_config) + "}$ & $\\num{" + str(worse_than_config_templated) + "}$ & $" + str(worst_quality) + "\degree$ & $" + str(worst_quality_templated) + "\degree$ & $"

        #if relative_worst_value > 1.3 or relative_worst_value < 0.7:
          #line += "\\color{red}"
        #line += str(relative_worst_value) + "$"

        #if worse_than_config_templated > worse_than_config:
          #in_use = True


      #if config_radius_edge_ratio != -1:
        #worse_than_config = float(get_value(filename, "Over config value", 1))
        #worse_than_config_templated = float(get_value(filename, "Over config value", 2))
        #worst_quality = round(float(get_value(filename, "max", 2)), 3)
        #worst_quality_templated = round(float(get_value(filename, "max", 4)), 3)

        #relative_worst_value = round(worst_quality_templated/config_radius_edge_ratio, 3)

        ##line += "$\\num{" + str(worse_than_config) + "}$ & $\\num{" + str(worse_than_config_templated) + "}$ & $" + str(worst_quality) + "$ & $" + str(worst_quality_templated) + "$ & $"

        ##if relative_worst_value > 1.3 or relative_worst_value < 0.7:
          ##line += "\\color{red}"

        ##line += str(relative_worst_value) + "$"

        #if worse_than_config_templated > worse_than_config:
          #in_use = True


      #line += " \\\\"

      #if in_use:
        #print(line)
        #print("\\hline")


check_quality("circle/quality")
