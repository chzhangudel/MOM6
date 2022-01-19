#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import copy as copy

def iou_py(diffu_py):
  diffu_py = diffu_py*2
#  diffu_py = np.float64(diffu_py)
#  print('diffu > \n',diffu_py[:,:,0])
#  print(diffu_py.dtype)
  return diffu_py

def iov_py(diffv_py):
  diffv_py = diffv_py*3
#  diffv_py = np.float64(diffv_py)
#  print('diffv > \n',diffv_py[:,:,0])
#  print(diffv_py.dtype)
  return diffv_py

