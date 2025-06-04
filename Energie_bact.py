import numpy as np
import matplotlib.pyplot as plt
import random as rd

def Etotal_glob(list_b,matrice_sucre):
    Etotal = 0
    eps= 1e-6  # petite valeur pour éviter la division par zéro
    sucre_total=np.sum(matrice_sucre) # total du sucre
    for bact in list_b:
        for i in range(matrice_sucre.shape[0]):
            for j in range(matrice_sucre.shape[1]):
                Rij_distance = ((bact.posx - j/matrice_sucre.shape[1])**2 + (bact.posy - i/matrice_sucre.shape[0])**2)**0.5
                if Rij_distance > 0 and sucre_total > 0  :  # éviter la division par zéro
                    Etotal += (matrice_sucre[i][j]/sucre_total)/ Rij_distance
    return Etotal

def gradglobal(list_b,matrice_sucre):
    epsilon = 1e-5  # petite valeur pour éviter la division par zéro
    for bact in list_b:
        Etot = Etotal_glob(list_b,matrice_sucre)
        bact.posx += epsilon
        Etot_prime = Etotal_glob(list_b,matrice_sucre)
        bact.posx -= epsilon
        grad_totx = (Etot_prime - Etot) / epsilon
        bact.posy += epsilon
        Etot_prime = Etotal_glob(list_b,matrice_sucre)
        grad_toty = (Etot_prime - Etot) / epsilon
        sum_grad = (grad_totx**2 + grad_toty**2)**0.5
        if sum_grad > 0:
            gradnumx = grad_totx / sum_grad
            gradnumy = grad_toty / sum_grad
            bact.gradnumx = gradnumx
            bact.gradnumy = gradnumy

def Etotal_local(bact,matrice_sucre):
    Etotal_local = 0
    sucre_total=np.sum(matrice_sucre) # total du sucre
    for i in range(matrice_sucre.shape[0]):
        for j in range(matrice_sucre.shape[1]):
            Etotal_local += (matrice_sucre[i][j]/sucre_total * 1)/ ((bact.posx - j/matrice_sucre.shape[1])**2 + (bact.posy - i/matrice_sucre.shape[0])**2)**0.5
    return Etotal_local

