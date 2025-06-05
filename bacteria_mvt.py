"""
Ce fichier décrit le mouvement des bactéries dans le bio-réacteur.
Les bactéries sont approximées par des points ayant des coordonnées x, y à tout moment t.
"""

import numpy as np
import random as rd

# ===================== PARAMÈTRES GLOBAUX =====================
VISION = 3                     # Distance de vision des bactéries ATTENTION METTRE QUE DES CHIFFRES IMPAIRS
LAMBDA_GRADIENT = 0.01          # Force du déplacement vers le gradient
LAMBDA_RANDOM = 0              # Intensité du mouvement brownien (aléatoire)
USE_GRADIENT = True             # Hypothèse 1 : déplacement selon gradient
USE_RANDOM = not USE_GRADIENT   # Hypothèse alternative

ATP_COST_MOVE = 10              # Coût d'ATP pour un déplacement
ATP_COST_MITOSIS = 300          # Coût d'ATP pour une division
ATP_THRESHOLD_MITOSIS = 200     # Seuil d'ATP pour se diviser

GLUCOSE_THRESHOLD_RESPIRATION = 0.1 # Seuil de glucose pour passer en respiration

MOLECULE_GLUC = 0.001
GLUCOSE_CONSUMPTION = 0.05      # Quantité de glucose consommée par itération
NB_MOLECULES_CONSOMMEES = GLUCOSE_CONSUMPTION/MOLECULE_GLUC
COEFF_CONSUMPTION_STATE = 30    # Coefficient pour la consommation de glucose en fonction de l'état
ATP_GAIN_EAT_FERMENTATION = 2 * NB_MOLECULES_CONSOMMEES  * COEFF_CONSUMPTION_STATE# Gain d'ATP en consommant du glucose
ATP_GAIN_EAT_RESPIRATION = 38 * NB_MOLECULES_CONSOMMEES  # Gain d'ATP en consommant du glucose en respiration
# ===================== CLASSE BACTERIA =====================
class Bacteria:
    def __init__(self, b_posx, b_posy):
        self.nom = 'E.coli'
        self.posx = b_posx
        self.posy = b_posy
        self.posmatx = 0
        self.posmaty = 0
        self.ATP = 100
        self.pop_threshold = 10        # Seuil de surpopulation
        self.death_pop_chance = 0.8    # Proba de mourir en cas de surpopulation
        self.death = False
        self.Etotal = 0
        self.gradx = 0
        self.grady = 0
        self.gradnumx = 0
        self.gradnumy = 0
        self.consommation_state = 'respiration'  # État de la bactérie
    # ===== Gradient perçu par la bactérie =====
    def Gradient(self, matrice):
        gluc_tot = np.sum(matrice)
        gradx = grady = 0

        for i in range(self.posmaty - VISION, self.posmaty + VISION):
            for j in range(self.posmatx - VISION, self.posmatx + VISION):
                if 0 <= i < matrice.shape[0] and 0 <= j < matrice.shape[1]:
                    dx = j - self.posmatx
                    dy = i - self.posmaty
                    dist = (dx**2 + dy**2)**0.5
                    if dist != 0:
                        coeff = matrice[i][j] / gluc_tot / dist**3
                        gradx += coeff * dx
                        grady += coeff * dy

        norm = (gradx**2 + grady**2)**0.5 or 1
        self.gradx = gradx / norm
        self.grady = grady / norm

    # ===== Mise à jour de la position =====
    def update_b_pos(self, matrice):
        taille = matrice.shape[0]
        self.posmatx = int(self.posx * taille)
        self.posmaty = int(self.posy * taille)

        self.Gradient(matrice)

        brownx = rd.uniform(-1, 1)
        browny = rd.uniform(-1, 1)

        newposx = self.posx + LAMBDA_RANDOM * brownx + LAMBDA_GRADIENT * self.gradx
        newposy = self.posy + LAMBDA_RANDOM * browny + LAMBDA_GRADIENT * self.grady

        if 0 <= newposx <= 1:
            self.posx = newposx
        if 0 <= newposy <= 1:
            self.posy = newposy

        self.ATP -= ATP_COST_MOVE

    # ===== Mort et mitose =====
    def update_death_and_mitosis(self, matrice, list_b):
        if self.ATP <= 0:
            self.death = True
            return

        # Compte les bactéries dans la même case
        nb_same_cell = sum(
            1 for bact in list_b if bact.posmatx == self.posmatx and bact.posmaty == self.posmaty
        )

        if nb_same_cell >= self.pop_threshold:
            if rd.choices([True, False], [self.death_pop_chance, 1 - self.death_pop_chance])[0]:
                self.death = True
                return

        # Mitose
        if nb_same_cell < self.pop_threshold and self.ATP > ATP_THRESHOLD_MITOSIS:
            dif = 0.001
            newx = self.posx + rd.uniform(-dif, dif)
            newy = self.posy + rd.uniform(-dif, dif)
            if 0 < newx < 1 and 0 < newy < 1:
                self.ATP -= ATP_COST_MITOSIS
                return Bacteria(newx, newy)

    # ===== Consommation de glucose =====
    def update_eat(self, matrice):
        if self.consommation_state == 'fermentation':    
            GLUCOSE_CONSUMPTION = COEFF_CONSUMPTION_STATE * GLUCOSE_CONSUMPTION
            if matrice[self.posmaty][self.posmatx] > 3 * GLUCOSE_CONSUMPTION:
                matrice[self.posmaty][self.posmatx] -= GLUCOSE_CONSUMPTION
                self.ATP += ATP_GAIN_EAT_FERMENTATION
        elif self.consommation_state == "respiration":
            if matrice[self.posmaty][self.posmatx] > 3 * GLUCOSE_CONSUMPTION:
                matrice[self.posmaty][self.posmatx] -= GLUCOSE_CONSUMPTION
                self.ATP += ATP_GAIN_EAT_RESPIRATION
        return matrice

    def update_state(self, matrice):
        moyenne_glucose = np.mean(matrice[self.posmaty - VISION:self.posmaty + VISION,
                                        self.posmatx - VISION:self.posmatx + VISION])
        if moyenne_glucose > GLUCOSE_THRESHOLD_RESPIRATION:
            self.consommation_state = 'fermentation'
        else:
            self.consommation_state = 'respiration'
