# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE FEM 1D
Title: config_packages.py
"""
import importlib
module_name = "data.data_ex1"
data = importlib.import_module(module_name)

# Basic packages:
import numpy as np 
import math
import matplotlib.pyplot as plt
from matplotlib import cm
import logging
