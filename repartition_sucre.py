"""
Ce script simule la diffusion de sucre dans une matrice 2D. Avec un flux modulable de glucose.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random

m_taille = 10 # Taille de la matrice
add_glucose = 10 # nb d'itérations avant l'ajout de glucose

def init_matrice(taille):
    """Initialise une matrice de taille donnée avec des zeroes."""
    return np.zeros((taille, taille))

#permet de mettre du sucre a un endroit donné de la matrice pour un certain amount et une certaine area
def sucre_input(mat,position,amount,area):
    x=position[0]
    y=position[1]
    for i in range(x-area,x+area):
        for j in range(y-area,y+area):
            if mat[i][j]+amount>1:
                mat[i][j]=1
            else:
                mat[i][j]+=amount
    return mat

#permet d'afficher la matrice avec une échelle de couleur grace à mathplotlib
def afficher_matrice(mat, titre=""):
    plt.imshow(mat, cmap='YlOrRd', vmin=0, vmax=1)
    plt.colorbar(label='Concentration de sucre')
    plt.title(titre)
    plt.axis('off')  # cache les axes
    plt.show()

#update la matrice de sucre en le diffusant dans les cases voisines en fonction du gradient de concentration
def update_sucre(old_mat):
    new_mat = np.copy(old_mat)
    D = 0.1  # coefficient de diffusion

    for i in range(m_taille):
        for j in range(m_taille):
            for di, dj in [(-1,0), (1,0), (0,-1), (0,1)]:
                if 0 <= i + di < m_taille and 0 <= j + dj < m_taille:
                    ni, nj = i + di, j + dj
                    gradient = old_mat[ni][nj] - old_mat[i][j]
                    flux = D * gradient      # Flick's LAw
                    new_mat[i][j] += flux
                    new_mat[ni][nj] -= flux  # Conservation de la masse
    # Clamping entre 0 et 1 pour éviter les valeurs absurdes
    new_mat = np.clip(new_mat, 0, 1)
    # Arrondie la matrice avec 4 décimale
    new_mat = np.round(new_mat,4)
    return new_mat


