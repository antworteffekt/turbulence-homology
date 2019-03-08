# -*- coding: utf-8 -*-
"""
Created on Tue May 29 15:45:23 2018

@author: licon
"""
import numpy as np
from homology.barcode import Barcode

class StableRankVectorizer(object):
    
    def __init__(self, embedding_dim=10, homology_dim=1,
                 persistence_type='additive', normalize=False,
                 threshold_min=0, threshold_max=5):
        """
        Vectorizer.
        """
#         self.input = input
        self.embedding_dim = embedding_dim
        self.homology_dim = homology_dim
        self.persistence_type = persistence_type
        self.normalize = normalize
        self.threshold_min = threshold_min
        self.threshold_max = threshold_max
    
    def get_params(self):

        return self.__dict__

    def transform(self, X):
        """
        Parameters
        ----------
        X: iterable containing barcodes to be vectorized.
        """
        self.evaluation_points = np.linspace(self.threshold_min, self.threshold_max,
                                             self.embedding_dim)
        X_vectorized = np.empty((len(X), self.embedding_dim))
        # Select stable rank function 
        if self.persistence_type == 'additive':
            stable_rank_func = getattr(Barcode, 'S_a')
        elif self.persistence_type == 'multiplicative':
            stable_rank_func = getattr(Barcode, 'S_p')
        else:
            raise ValueError("Persistence type not defined.")
        
        if self.persistence_type == 'multiplicative' and self.homology_dim < 1:
            raise ValueError("Multiplicative persistence not defined for dim < 1.")
        
        i = 0
        for barcode in X:
            X_vectorized[i] = stable_rank_func(barcode, self.evaluation_points, 
                                               self.homology_dim)
            i += 1
        return X_vectorized

