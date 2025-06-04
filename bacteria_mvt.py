"""
Ce fichier décrit le mouvement des bactéries dans le bio-réacteur
Les bactéries sont approximés selon des points ayant des coordonées x,y à tout moment t
"""
import numpy as np
import random as rd
#from main import *
vision_b = 5 #vision des bactéries
lambda1 = 0.01 #Vitesse de déplacement des bactéries selon le gradient
lambda2 = 0.0025 #Amplitude de l'aléatoire dans le déplacement des bactéries
hypothese1 = True#les bacteries se déplacent selon le gradient de concentration de glucose
hypothese2 = not hypothese1 #quantité de glucose consommée par les bactéries à chaque itération
kcoutatp = 20 #coût d'ATP pour se déplacer
kgainatp = 38 #gain d'ATP pour manger du glucose
kcoutmitose = 120 #coût d'ATP pour se diviser
kthresholdmitose = 200 #seuil d'ATP pour se diviser
kglucose = 0.1 #quantité de glucose consommée par les bactéries à chaque itération

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
        self.Etotal = 0
        self.gradx = 0
        self.grady = 0
    def Gradient(self,matrice):
        gluc_tot = np.sum(matrice)
        gradx=0
        grady=0  
        for i in range(self.posmaty-vision_b, self.posmaty+vision_b):
            for j in range(self.posmatx-vision_b, self.posmatx+vision_b):
                if 0<i<matrice.shape[0] and 0<j<matrice.shape[1]:
                    hypothenuse=((self.posmatx - j)**2 + (self.posmaty - i)**2)**0.5
                    if hypothenuse != 0:  # éviter la division par zéro
                        delta_distance_x= (j - self.posmatx) / hypothenuse**3
                        delta_distance_y= (i - self.posmaty) / hypothenuse**3
                        gradx += matrice[i][j]/gluc_tot * delta_distance_x
                        grady += matrice[i][j]/gluc_tot * delta_distance_y
        gradnx = gradx/(gradx**2+grady**2)**0.5
        gradny = grady/(gradx**2+grady**2)**0.5
        self.gradx = gradnx
        self.grady = gradny

    def update_b_pos(self, matrice):
        m_taille = matrice.shape[0]
        self.posmatx = int(m_taille * self.posx)
        self.posmaty = int(m_taille * self.posy)
        self.Gradient(matrice)
        if hypothese1:
            newposx = self.posx + lambda1 * self.gradx
            newposy = self.posy + lambda1 * self.grady
        #hypothese2 : les bacteries se déplacent avec plus de chance vers le glucose mais pas forcement
        if hypothese2:
        #marche pas pour l'instant à réfléchir comment faire pour que ca marche
        # composante brownienne aléatoire
            brownx = rd.uniform(-1, 1)
            browny = rd.uniform(-1, 1)
        # ajout d’un biais vers le gradient
            newposx = self.posx + lambda2 * brownx + lambda1 * self.gradx
            newposy = self.posy + lambda2 * browny + lambda1 * self.grady
        #debug des murs : si la posx ou posy sort de la matrice on n'avance pas
        if 0 <= newposx <= 1 : self.posx=newposx
        if 0 <= newposy <= 1 : self.posy=newposy
        #update de la position de la bactérie dans la matrice
        self.ATP -= kcoutatp
        
    def update_death_and_mitosis(self,matrice,list_b):
        #mort de la bactérie si il y a trop de bactéries sur la case
        nb_bacteries_case = 0
        for bact in list_b:
            if bact.posmatx == self.posmatx and bact.posmaty == self.posmaty:
                nb_bacteries_case += 1
        if nb_bacteries_case >= self.death_threshold:
            if rd.choices([True, False], [self.death_chance, 1-self.death_chance]):
                self.death = True
        #mitosis si il y a assez d'ATP et pas trop de bactéries sur la case
        if nb_bacteries_case < self.death_threshold and self.ATP > kthresholdmitose:    
                # On crée une nouvelle bactérie à une position légèrement différente
                new_bacteria = Bacteria(self.posx,self.posy)
                self.ATP -= kcoutatp
                return new_bacteria
            
    def update_eat(self, matrice):
        # On mange le glucose de la case de la matrice
        if matrice[self.posmaty][self.posmatx] > kglucose*3:   # On mange que si il y a assez de glucose    
            matrice[self.posmaty][self.posmatx] -= kglucose # Consommation de glucose 
            self.ATP += kgainatp # On gagne de l'ATP
        return matrice 

        
