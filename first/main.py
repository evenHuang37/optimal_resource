#!/usr/bin/env python3
# -*- coding:utf-8 -*-
__author__ = 'Even Huang'

import numpy
import matplotlib.pyplot as plt
from first.generate_ER import *
from first.SIS import *
import re



if __name__ == "__main__":
    # degree_amount = {}
    # degree_amount[degree] = degree_amount.get(degree,0) + 1
    # sorted = sorted(degree_amount.items(),key = lambda a:a[0])

    beta = 0.1
    miu = 0.2
    t = 50
    repeat_times = 10
    alpha = 0
    p_0 = 0.05
    N = 100
    edge_p = 0.06

    sis = SIS_model(N, edge_p, beta, miu, t, p_0, "ER_network", "density", repeat_times)
    sis.Graw_SIS()