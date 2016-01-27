# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 10:32:00 2014

@author: Kyle Ellefsen
"""
from collections import OrderedDict


name='Detect Puffs'
base_dir='detect_puffs'
date='11/15/2015'
dependencies = ['skimage', 'leastsqbound', 'matplotlib','PyOpenGL']

menu_layout = {'Detect Puffs':  OrderedDict([    ('Threshold Cluster',   ['threshold_cluster', 'threshold_cluster.gui']),
				                          ('Simulate Puffs'   ,   ['puff_simulator.puff_simulator', 'simulate_puffs.gui']),
                                                 ('Load .flika file' ,   ['threshold_cluster', 'load_flika_file_gui']),
                                                 ('Docs'             ,   ['threshold_cluster', 'launch_docs'])
				                    ])
               }