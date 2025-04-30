"""
Ce fichier décrit le mouvement des bactéries dans le bio-réacteur
Les bactéries sont approximés selon des points ayant des coordonées x,y à tout moment t
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random
lambdaa = 0.01
class Bacteria:
    def __init__(self, b_posx,b_posy):
        self.nom = 'E.coli'
        self.posx = b_posx
        self.posy = b_posy
        self.death_rate = 0
        self.mitose_rate = 0
        self.ATP_bank = 0
    def update_b_pos(self, matrice):
        index_casex = int(matrice.shape[0] * self.posx)
        index_casey = int(matrice.shape[1] * self.posy)
        ratio_casex = matrice.shape[0] * self.posx - index_casex
        ratio_casey = matrice.shape[1] * self.posy - index_casey
        gradientx = matrice[index_casex+1][index_casey] - matrice[index_casex-1][index_casey]
        gradienty = matrice[index_casex][index_casey+1] - matrice[index_casex][index_casey-1]
        self.posx = self.posx + (lambdaa * gradientx)
        self.posy = self.posy + (lambdaa * gradienty)

"""
list_b = []
for x in list_b:
    x.update_b_pos(matrice)
"""