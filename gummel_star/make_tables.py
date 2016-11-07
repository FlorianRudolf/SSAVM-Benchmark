#!/usr/bin/python

import sys
from dolfin import *
from math import pi
import os


def ftos(x):
  s = "%.1E" % x
  [a,b] = s.split("E")
  if b[1] == '0':
    b = "-" + b[2] + "\phantom{0}"
    #b[1] = b[0]
    #b[0] = '0'
    #b = b.replace("0", "\phantom{0}")
  return a + " \\times 10^{" + b + "}"


fn_dict = { "results_nonsym_gs_c0.016384" : "Non-symmetric & $109$",
            "results_reflrotsym_gs_c0.016384" : "Reflective-rotationally symmetric & $130$",
            "results_rotsym_gs_c0.016384" : "Rotationally symmetric & $120$",

            "results_nonsym_gs_c0.001024" : "Non-symmetric & $1758$",
            "results_reflrotsym_gs_c0.001024" : "Reflective-rotationally symmetric & $1750$",
            "results_rotsym_gs_c0.001024" : "Rotationally symmetric & $1740$",

            "results_nonsym_gs_c0.000128" : "Non-symmetric & $13921$",
            "results_reflrotsym_gs_c0.000128" : "Reflective-rotationally symmetric & $14160$",
            "results_rotsym_gs_c0.000128" : "Rotationally symmetric & $13850$",

            "results_nonsym_gs_c1.6e-05" : "Non-symmetric & $111244$",
            "results_reflrotsym_gs_c1.6e-05" : "Reflective-rotationally symmetric & $111660$",
            "results_rotsym_gs_c1.6e-05" : "Rotationally symmetric & $111260$",

           }


filenames=["results_nonsym_gs_c0.016384", "results_rotsym_gs_c0.016384", "results_reflrotsym_gs_c0.016384",
           "results_nonsym_gs_c0.001024", "results_rotsym_gs_c0.001024", "results_reflrotsym_gs_c0.001024",
           "results_nonsym_gs_c0.000128", "results_rotsym_gs_c0.000128", "results_reflrotsym_gs_c0.000128",
           "results_nonsym_gs_c1.6e-05", "results_rotsym_gs_c1.6e-05", "results_reflrotsym_gs_c1.6e-05"]

dual_table = [[0 for x in range(6)] for y in range(len(filenames))]
cn = 0

for filename in filenames:
  file = open(filename, "r")

  table  = [[0 for x in range(5)] for y in range(5)]
  dual_table[cn][0] = filename

  for line in file:
    line = line[:-1]
    tmp = line.split(";")

    if tmp[1] != "dual" and tmp[3] != "dual":
      table[int(tmp[0])][int(tmp[2])] = float(tmp[4])

    #print(tmp)
    if (tmp[0] == tmp[2]) and (tmp[1] == "normal") and (tmp[3] == "dual"):
      dual_table[cn][int(tmp[0])+1] = float(tmp[4])

  #print table

  of = open("tables/" + filename + "_normal", "w")

  #of.write("\\begin{tabular}{|l|c|c|c|c|c|c|} \n")
  #of.write("  \\hline\n")
  #of.write("  & \\multilinevertical{No \\\\ rotation} & \\multilinevertical{$72\\degree$ \\\\ Rotation} & \multilinevertical{$144\\degree$ \\\\ Rotation} & \\multilinevertical{$216\\degree$ \\\\ Rotation} & \\multilinevertical{$288\\degree$ \\\\ Rotation$\\phantom{0}$} \\\\ \n")
  #of.write("  \\hline \n")

  for r in range(5):

    if r == 0:
      line = "No rotation"
    else:
      line = "Rotation $" + str(72*r) + "\\degree$"

    #line += " & "

    for c in range(5):
      line += " & "

      if (r == c):
        line += "$\color{OliveGreen} 0$"
      elif (r < c):
        line += "$"
        if (abs(table[r][c]) < 1e-8):
          line += "\color{OliveGreen} "
        if (abs(table[r][c]) > 1e-2):
          line += "\color{red} "
        line += ftos(table[r][c]) + "$"

    line += "\\\\ \n"
    of.write(line)
    of.write("  \\hline\n")


  cn = cn+1

of = open("tables/dual", "w")
for entry in dual_table:
  line = fn_dict[entry[0]]
  for c in range(1):
    line += " & $"
    if (abs(entry[c+1]) < 1e-8):
      line += "\color{OliveGreen} "
    if (abs(entry[c+1]) > 1e-2):
      line += "\color{red} "
    line += ftos(entry[c+1]) + "$"

  line += "\\\\ \n"

  of.write(line)
  of.write("\\hline\n")


  #of.write("\end{tabular} \n")


    #key = tmp[0:4];
    #distance = tmp[4]

    #values[key] = distance

    #print tmp
    #print key