"""
Ce fichier décrit le mouvement des bactéries dans le bio-réacteur
Les bactéries sont approximés selon des points ayant des coordonées x,y à tout moment t
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random

class Bacteria:
    def __init__(self, b_posx,b_posy):
        self.nom = 'E.coli'
        self.posx = b_posx
        self.posy = b_posy
        self.death_rate = 0
        self.mitose_rate = 0
        self.ATP_bank = 0
    def update_b_pos(self, matrice):
        


list_b = []
for x in list_b:
    x.update_b_pos(matrice)