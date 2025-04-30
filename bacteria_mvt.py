"""
Ce fichier décrit le mouvement des bactéries dans le bio-réacteur
Les bactéries sont approximés selon des points ayant des coordonées x,y à tout moment t
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random
lambdaa = 0.1
class Bacteria:
    def __init__(self, b_posx,b_posy):
        self.nom = 'E.coli'
        self.posx = b_posx
        self.posy = b_posy
        self.death_rate = 0
        self.mitose_rate = 0
        self.ATP_bank = 0
    def update_b_pos(self, matrice):
        m_taille = matrice.shape[0]
        index_casex = int(m_taille * self.posx)
        index_casey = int(m_taille * self.posy)
        ratio_casex = m_taille * self.posx - index_casex
        ratio_casey = m_taille * self.posy - index_casey
        if 0 <= index_casex-1 < m_taille and 0 <= index_casex+1 < m_taille:
            gradientx = matrice[index_casex+1][index_casey] - matrice[index_casex-1][index_casey]
            self.posx = self.posx + (lambdaa * ratio_casex * gradientx)
        if 0 <= index_casey-1 < m_taille and 0 <= index_casey+1 < m_taille:
            gradienty = matrice[index_casex][index_casey+1] - matrice[index_casex][index_casey-1]
            self.posy = self.posy + (lambdaa * ratio_casey * gradienty)

"""
list_b = []
for x in list_b:
    x.update_b_pos(matrice)
"""