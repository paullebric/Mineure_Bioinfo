"""
Ce fichier décrit le mouvement des bactéries dans le bio-réacteur
Les bactéries sont approximés selon des points ayant des coordonées x,y à tout moment t
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random
lambda1 = 0.001
lambda2 = 0.00001
hyposthese1 = True  #les bacteries se déplacent selon le gradient de concentration de glucose
hyposthese2 = False #les bacteries se déplacent proportionnellement vers le glucose mais pas forcement
class Bacteria:
    def __init__(self, b_posx,b_posy):
        self.nom = 'E.coli'
        self.posx = b_posx
        self.posy = b_posy
        self.posmaty = 0
        self.posmatx = 0
        self.ATP = 100
        self.death_threshold = 10
        self.death_chance = 0.5
        self.death = False
    def update_b_pos(self, matrice):
        m_taille = matrice.shape[0]
        index_casex = int(m_taille * self.posx)
        index_casey = int(m_taille * self.posy)
        # ratio_casex = m_taille * self.posx - index_casex
        # ratio_casey = m_taille * self.posy - index_casey
        index_casedroite = index_casex + 1
        index_casegauche = index_casex - 1
        index_casehaut = index_casey + 1
        index_casebas = index_casey - 1
        if index_casegauche <0: index_casex = m_taille-1
        if index_casebas <0: index_casey = m_taille-1
        if index_casehaut >= m_taille: index_casehaut = 0
        if index_casedroite >= m_taille: index_casedroite = 0            
        gradientx = matrice[index_casedroite][index_casey] - matrice[index_casegauche][index_casey]
        gradienty = matrice[index_casex][index_casehaut] - matrice[index_casex][index_casebas]
        gradn = np.sqrt(gradientx**2 + gradienty**2)
        if gradn == 0:
            gradn = 1e-10
        delta_gradientx = gradientx/gradn
        delta_gradienty = gradienty/gradn
        if hyposthese1:
            self.posx = self.posx + lambda1 * delta_gradientx+lambda2 * random.uniform(-1, 1)
            self.posy = self.posy + lambda1 * delta_gradienty+lambda2 * random.uniform(-1, 1)

        self.posmatx = int(m_taille * self.posx)
        self.posmaty = int(m_taille * self.posy)
    def update_death_and_mitosis(self,matrice,list_b):
        nb_bacteries_case = 0
        for bact in list_b:
            if bact.posmatx == self.posmatx and bact.posmaty == self.posmaty:
                nb_bacteries_case += 1
        if nb_bacteries_case >= self.death_threshold:
            if random.choices([True, False], [self.death_chance, 1-self.death_chance]):
                self.death = True
                
        
        #cette fonction permet de manger du glucose
        

"""
list_b = []
for x in list_b:
    x.update_b_pos(matrice)
"""