# -*- coding: utf-8 -*-
"""
Created on Tue May  5 13:13:54 2020

@author: OleK
"""

import numpy as n

Re=6371.0
latdeg2km=n.sin(n.pi/180.0)*Re
londeg2km=n.pi*Re*n.cos(n.pi*69.0/180.0)/180.0#65.122785 
