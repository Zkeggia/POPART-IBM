#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 12:20:53 2021

@author: fra
"""

from treetime import TreeTime,plot_vs_years
from treetime.utils import parse_dates

dates=parse_dates('Outputs/output_sample_times.tsv', name_col='#label', date_col='date')
tt=TreeTime(tree='Outputs/sequences.PHYLIP.treefile', aln='Outputs/sequences.PHYLIP',  dates=dates)
tt.run(root='best')

plot_vs_years(tt)