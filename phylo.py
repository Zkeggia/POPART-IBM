#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 13:07:49 2021

@author: fra
"""

import pandas as pd, numpy as np
import argparse
from Bio import Phylo

if __name__=="__main__":
    tree = Phylo.read('Outputs/phylogenetic.newick','newick')
    
    Phylo.draw(tree)
    
