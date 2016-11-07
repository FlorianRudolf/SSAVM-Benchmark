#!/usr/bin/python

import sys
from dolfin import *
from math import pi
import os

binary = "../../tools/build_release/compare_solutions"
path = "results"
filenames = ["gs_c0.016384", "gs_c0.001024", "gs_c0.000128", "gs_c1.6e-05"]
#filenames = ["gs_c0.016384"]
#filenames = ["gs_c0.001024"]
symmetries = ["nonsym", "rotsym", "reflrotsym"]
#symmetries = ["nonsym", "reflrotsym"]

for fn in filenames:
  for sym in symmetries:
    f = open("results_" + sym + "_" + fn, 'w')

    for rot1 in [0,1,2,3,4,5,6,7,8,9]:
      for rot2 in range(rot1, 10):

        dual1 = False
        dual2 = False

        if rot1 >= 5:
          r1 = rot1-5
          dual1 = True
        else:
          r1 = rot1

        if rot2 >= 5:
          r2 = rot2-5
          dual2 = True
        else:
          r2 = rot2

        filename1 = path + "/solution_" + ("dual_" if dual1 else "") + str(r1) + "_" + fn + "_" + sym + ".pvd"
        filename2 = path + "/solution_" + ("dual_" if dual2 else "") + str(r2) + "_" + fn + "_" + sym + ".pvd"

        result = float(os.popen( binary + " " + filename1 + " " + filename2 + " --rotation1 " + str(r1) + " --rotation2 " + str(r2) + (" --dual1" if dual1 else "") + (" --dual2" if dual2 else "")).read())

        string = str(r1) + ";" + ("dual" if dual1 else "normal") + ";" + str(r2) + ";" + ("dual" if dual2 else "normal") + ";" + str(result)
        print(fn + ";" + sym + ";" + string)
        f.write(string + "\n")

    print("\n\n\n")


#for fn in filenames:
  #for sym in symmetries:
    #for rot1 in [0,1,2,3,4]:
      #for rot2 in range(rot1, 5):
        #filename1 = path + "/solution_dual_" + str(rot1) + "_" + fn + "_" + sym + ".pvd"
        #filename2 = path + "/solution_dual_" + str(rot2) + "_" + fn + "_" + sym + ".pvd"

        #result = float(os.popen( binary + " " + filename1 + " " + filename2 + " --rotation1 " + str(rot1) + " --rotation2 " + str(rot2) ).read())
        #print("DUAL " + fn + " " + sym + " " + str(rot1) + " " + str(rot2) + "    " + str(result))


 #results/solution_0_gs_c0.016384_nonsym.pvd results/solution_dual_1_gs_c0.016384_nonsym.pvd --rotation1 0 --rotation2 1 --dual2