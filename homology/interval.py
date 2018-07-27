# -*- coding: utf-8 -*-
"""
Created on Tue May 29 15:44:52 2018

@author: licon
"""

from __future__ import print_function

class Interval:
    
    def __init__(self, start, end=None, is_finite=True):
        self.start = start
        self.end = end
        self.is_finite = is_finite
    
    
    def __repr__(self):
        return(
            'Persistence interval:\n\tstart={}, end={}'
            .format(
                self.start,
                self.end
            ))
        
        
