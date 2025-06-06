"""
Ce fichier décrit le mouvement des bactéries dans le bio-réacteur.
Les bactéries sont approximées par des points ayant des coordonnées x, y à tout moment t.
"""

import numpy as np
import random as rd

# ===================== PARAMÈTRES GLOBAUX =====================
VISION = 3     # Distance de vision des bactéries ATTENTION METTRE QUE DES CHIFFRES IMPAIRS
LAMBDA_GRADIENT = 0.07          # Force du déplacement vers le gradient
LAMBDA_RANDOM = 0.001          # Intensité du mouvement brownien (aléatoire)

NOMBRE_LIMITE_MITOSE = 5  # Nombre limite de mitoses pour une bactérie

ATP_COST_MOVE = 10              # Coût d'ATP pour un déplacement
ATP_COST_MITOSIS = 1000          # Coût d'ATP pour une division
ATP_THRESHOLD_MITOSIS = 1200    # Seuil d'ATP pour se diviser
MITOSE_CHANCE = 0.25 # Proba de se diviser si conditions remplies

GLUCOSE_THRESHOLD_FERMENTATION = 0.15 # Seuil de glucose pour passer en respiration

MOLECULE_GLUC = 0.001
GLUCOSE_CONSUMPTION = 0.00125     # Quantité de glucose consommée par itération
NB_MOLECULES_CONSOMMEES = GLUCOSE_CONSUMPTION/MOLECULE_GLUC
COEFF_CONSUMPTION_STATE = 30    # Coefficient pour la consommation de glucose en fonction de l'état
ATP_GAIN_EAT_FERMENTATION = 2 * NB_MOLECULES_CONSOMMEES  * COEFF_CONSUMPTION_STATE# Gain d'ATP en consommant du glucose
ATP_GAIN_EAT_RESPIRATION = 38 * NB_MOLECULES_CONSOMMEES  # Gain d'ATP en consommant du glucose en respiration

ferment_tox = 0.001
respi_tox = 0.002
MAX_TOXINE = 0.8  # Concentration maximale de toxine dans une case
# ===================== CLASSE BACTERIA =====================
class Bacteria:
    def __init__(self, b_posx, b_posy):
        self.nom = 'E.coli'
        self.posx = b_posx
        self.posy = b_posy
        self.posmatx = 0
        self.posmaty = 0
        self.ATP = 100
        self.pop_threshold = 5        # Seuil de surpopulation
        self.death_pop_chance = 0.8    # Proba de mourir en cas de surpopulation
        self.death = False
        self.Etotal = 0
        self.gradx = 0
        self.grady = 0
        self.gradnumx = 0
        self.gradnumy = 0
        self.consommation_state = 'respiration'  # État de la bactérie
        self.mitose_count = 0  # Compteur de mitoses
        self.Rij3 = 0  # Variable pour la mise à jour de la position
        self.age = 100
    # ===== Gradient perçu par la bactérie =====
    def Gradient(self, matrice_sucre):
        gluc_tot = np.sum(matrice_sucre)
        gradx = grady = 0

        for i in range(self.posmaty - VISION, self.posmaty + VISION+1):
            for j in range(self.posmatx - VISION, self.posmatx + VISION+1):
                if 0 <= i < matrice_sucre.shape[0] and 0 <= j < matrice_sucre.shape[1]:
                    dx = (j+0.5)/ matrice_sucre.shape[1] - self.posx
                    dy = (i+0.5)/ matrice_sucre.shape[0] - self.posy
                    dist = (dx**2 + dy**2)**0.5
                    if dist != 0:
                        coeff = (matrice_sucre[i][j] / gluc_tot) / dist**3
                        gradx += coeff * dx
                        grady += coeff * dy


        norm = (gradx**2 + grady**2)**0.5
        self.gradx = gradx / norm
        self.grady = grady / norm
    def update_Rij3(self, matrice_sucre):
        for i in range(self.posmaty - VISION, self.posmaty + VISION+1):
            for j in range(self.posmatx - VISION, self.posmatx + VISION+1):
                if 0 <= i < matrice_sucre.shape[0] and 0 <= j < matrice_sucre.shape[1]:
                    dx = (j+0.5)/ matrice_sucre.shape[1] - self.posx
                    dy = (i+0.5)/ matrice_sucre.shape[0] - self.posy
                    dist = (dx**2 + dy**2)**0.5
                    self.Rij3 = 1/(dist**3)

    def Gradient_concerte(self, matrice_sucre, list_b):
        gluc_tot = np.sum(matrice_sucre)
        gradx = grady = 0
        sum_Rij3 = 0
        sum_Sxj2 = 0
        sum_Syj2 = 0
        for i in range(self.posmaty - VISION, self.posmaty + VISION + 1):
                    for j in range(self.posmatx - VISION, self.posmatx + VISION + 1):
                        if 0 <= i < matrice_sucre.shape[0] and 0 <= j < matrice_sucre.shape[1]:
                            gluc_pourcentage = matrice_sucre[i][j] / gluc_tot
                            sum_Sxj2 -= gluc_pourcentage*((j + 0.5)/matrice_sucre.shape[0]-self.posx)
                            sum_Syj2 -= gluc_pourcentage*((i + 0.5)/matrice_sucre.shape[1]-self.posy)

        for bact in list_b:
                                if self.posmaty - VISION < bact.posmaty < self.posmaty + VISION + 1 and self.posmatx - VISION < bact.posmatx < self.posmatx + VISION + 1:
                                    sum_Rij3 += bact.Rij3
        gradx = sum_Sxj2*sum_Rij3
        grady = sum_Syj2*sum_Rij3
        norm = (gradx**2 + grady**2)**0.5
        self.gradx = gradx / norm
        self.grady = grady / norm
        #print("Gradient concerté calculé")
    # ===== Mise à jour de la position =====
    def update_b_pos(self, matrice_sucre, CONCERTE, list_b):
        taille = matrice_sucre.shape[0]
        self.posmatx = int(self.posx * taille)
        self.posmaty = int(self.posy * taille)
        if CONCERTE :
            self.Gradient_concerte(matrice_sucre, list_b)
            #print("Gradient concerté")
        else :
            self.Gradient(matrice_sucre)

        brownx = rd.uniform(-1, 1)
        browny = rd.uniform(-1, 1)

        newposx = self.posx + LAMBDA_RANDOM * brownx + LAMBDA_GRADIENT * self.gradx
        newposy = self.posy + LAMBDA_RANDOM * browny + LAMBDA_GRADIENT * self.grady
        #print("Nouvelle position :", newposx, newposy)
        if 0 <= newposx <= 1:
            self.posx = newposx
        if 0 <= newposy <= 1:
            self.posy = newposy

        self.ATP -= ATP_COST_MOVE

        
    # ===== Mort et mitose =====
    def update_death_and_mitosis(self, mat_tox, list_b):
        self.age -= 1
        if self.age <= 0:
            self.death = True
            #print("bactérie morte (âge)")
            return
        if self.ATP <= 0 or self.mitose_count >= NOMBRE_LIMITE_MITOSE:
            self.death = True
            #print("bactérie morte")
            #print(self.ATP)
            #print(self.mitose_count)
            return
        if mat_tox[self.posmaty][self.posmatx] > MAX_TOXINE:
            self.death = True
            #print("bactérie morte (toxine)")
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
            if rd.random() < MITOSE_CHANCE:
                dif = 0.001
                newx = self.posx + rd.uniform(-dif, dif)
                newy = self.posy + rd.uniform(-dif, dif)
                if 0 < newx < 1 and 0 < newy < 1:
                    self.ATP -= ATP_COST_MITOSIS
                    self.mitose_count += 1
                    new_bact = Bacteria(newx, newy)
                    new_bact.consommation_state = self.consommation_state
                    
                    return new_bact

    # ===== Consommation de glucose =====
    def update_eat(self, matrice_sucre, mat_tox, GLUCOSE_CONSUMPTION=GLUCOSE_CONSUMPTION):
        if self.consommation_state == 'fermentation':    
            ferment_GLUCOSE_CONSUMPTION = COEFF_CONSUMPTION_STATE * GLUCOSE_CONSUMPTION
            if matrice_sucre[self.posmaty][self.posmatx] > 3 * ferment_GLUCOSE_CONSUMPTION:
                matrice_sucre[self.posmaty][self.posmatx] -= ferment_GLUCOSE_CONSUMPTION
                self.ATP += ATP_GAIN_EAT_FERMENTATION
                mat_tox[self.posmaty][self.posmatx] += ferment_tox  # Ajoute une molécule de glucose dans la matrice toxine
        elif self.consommation_state == "respiration":
            if matrice_sucre[self.posmaty][self.posmatx] > 3 * GLUCOSE_CONSUMPTION:
                matrice_sucre[self.posmaty][self.posmatx] -= GLUCOSE_CONSUMPTION
                self.ATP += ATP_GAIN_EAT_RESPIRATION
                mat_tox[self.posmaty][self.posmatx] += respi_tox
        return [matrice_sucre, mat_tox]

    def update_state(self, matrice_sucre):
        n_case = 0
        sum_gluc = 0
        for a in range(self.posmaty - VISION, self.posmaty + VISION+1):
            for b in range(self.posmatx - VISION, self.posmatx + VISION+1):
                if 0 <= a < matrice_sucre.shape[0] and 0 <= b < matrice_sucre.shape[1]:
                        sum_gluc+= matrice_sucre[a][b]
                        n_case+=1
        moyenne_glucose = sum_gluc / n_case if n_case > 0 else 0
        if moyenne_glucose > GLUCOSE_THRESHOLD_FERMENTATION:
            self.consommation_state = 'fermentation'
        else:
            self.consommation_state = 'respiration'
            
