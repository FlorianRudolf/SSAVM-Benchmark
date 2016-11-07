#!/usr/bin/python

import os


def file_valid( filename ):
  if os.stat(filename).st_size == 0:
    return False

  return True


def get_value( filename, value, number = 1 ):
  f = open(filename, 'r')
  i = 1
  for line in f:
    if line.find(value) != 0: continue
    if i == number:
      tokens = line.split()
      return tokens[-1]
    else:
      i = i+1

