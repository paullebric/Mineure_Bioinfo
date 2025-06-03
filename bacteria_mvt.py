"""
Ce fichier décrit le mouvement des bactéries dans le bio-réacteur
Les bactéries sont approximés selon des points ayant des coordonées x,y à tout moment t
"""

import numpy as np
import random as rd
#from main import *
lambda1 = 0.5 #amplitude de la composante brownienne aléatoire
lambda2 = 0.5 #vitesse de déplacement des bactéries
hypothese1 = True  #les bacteries se déplacent selon le gradient de concentration de glucose
hypothese2 = False #quantité de glucose consommée par les bactéries à chaque itération
kcoutatp = 20 #coût d'ATP pour se déplacer
kgainatp = 38 #gain d'ATP pour manger du glucose
kcoutmitose = 120 #coût d'ATP pour se diviser
kthresholdmitose = 200 #seuil d'ATP pour se diviser

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
    def Etotal_Gradient(self,matrice):
        sum_Elocal = 0
        sum_gradx = 0
        sum_grady = 0
        for a in range(matrice.shape[0]):
            for b in range(matrice.shape[1]):
                sum_Elocal += matrice[a][b]/(((a-self.posmaty)**2 + (b-self.posmatx)**2)**0.5)
                sum_gradx += (b-self.posmatx)*matrice[a][b]/(((a-self.posmaty)**2 + (b-self.posmatx)**2)**1.5)
                sum_grady += (a-self.posmaty)*matrice[a][b]/(((a-self.posmaty)**2 + (b-self.posmatx)**2)**1.5)
        print("sum", sum_gradx, sum_grady)
        self.Etotal = matrice[self.posmaty,self.posmatx]*sum_Elocal
        self.gradx = sum_gradx
        self.grady = sum_grady
        
    def update_b_pos(self, matrice):
        m_taille = matrice.shape[0]
        # index_casex = int(m_taille * self.posx)
        # index_casey = int(m_taille * self.posy)
        # # ratio_casex = m_taille * self.posx - index_casex
        # # ratio_casey = m_taille * self.posy - index_casey
        # index_casedroite = index_casex + 1
        # index_casegauche = index_casex - 1
        # index_casehaut = index_casey + 1
        # index_casebas = index_casey - 1
        # if index_casegauche <0: index_casex +=1
        # if index_casebas <0: index_casey +=1
        # if index_casehaut >= m_taille: index_casehaut -=1 
        # if index_casedroite >= m_taille: index_casedroite -=1            
        # gradientx = matrice[index_casey][index_casedroite] - matrice[index_casey][index_casegauche]
        # gradienty = matrice[index_casehaut][index_casex] - matrice[index_casebas][index_casex]
        # gradn = np.sqrt(gradientx**2 + gradienty**2)
        # if gradn == 0:
        #     gradn = 1e-10
        # delta_gradientx = gradientx/gradn
        # delta_gradienty = gradienty/gradn
        #hypothese1 : les bacteries se déplacent selon le gradient de concentration de glucose plus au moins vite aléatoirement
        self.Etotal_Gradient(matrice)
        if hypothese1:
            print("grad",self.gradx,self.grady)
            newposx = self.posx + lambda1 * self.gradx+lambda2
            newposy = self.posy + lambda1 * self.grady+lambda2
            print(newposx,newposy)
        #hypothese2 : les bacteries se déplacent avec plus de chance vers le glucose mais pas forcement
        if hypothese2:
        #marche pas pour l'instant à réfléchir comment faire pour que ca marche
        # composante brownienne aléatoire
            brownx = rd.uniform(-1, 1)
            browny = rd.uniform(-1, 1)
        # ajout d’un biais vers le gradient
            newposx = self.posx + lambda1 * brownx + lambda2 * self.gradx
            newposy = self.posy + lambda1 * browny + lambda2 * self.grady
        #debug des murs : si la posx ou posy sort de la matrice on n'avance pas
        if 0 <= newposx <= 1 : self.posx=newposx
        if 0 <= newposy <= 1 : self.posy=newposy
        #update de la position de la bactérie dans la matrice
        self.posmatx = int(m_taille * self.posx)
        self.posmaty = int(m_taille * self.posy)
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

        
