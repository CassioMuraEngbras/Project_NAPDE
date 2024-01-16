# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: NAPDE
Title: config_packages.py
"""
import importlib
module_name = "Data.data_ex3"
data = importlib.import_module(module_name)

# Basic packages:
import numpy as np 
import math
import matplotlib.pyplot as plt
from matplotlib import cm
