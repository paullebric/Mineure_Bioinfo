"""
Ce script simule la diffusion de sucre dans une matrice_sucre 2D. Avec un flux modulable de glucose.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random


def init_matrice_sucre(taille,sucre_base):
    #Initialise une matrice_sucre de taille donnée avec des zeroes.
    return np.full((taille, taille), float(sucre_base))


#permet de mettre du sucre a un endroit donné de la matrice_sucre pour un certain amount et une certaine area
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

#permet d'afficher la matrice_sucre avec une échelle de couleur grace à mathplotlib
def afficher_matrice_sucre(mat, titre=""):
    plt.imshow(mat, cmap='YlOrRd', vmin=0, vmax=1)
    plt.colorbar(label='Concentration de sucre')
    plt.title(titre)
    plt.axis('off')  # cache les axes
    plt.show()

#update la matrice_sucre de sucre en le diffusant dans les cases voisines en fonction du gradient de concentration
def update_sucre(old_mat):
    new_mat = np.copy(old_mat)
    D = 0.1  # coefficient de diffusion
    m_taille = old_mat.shape[0]
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
    # Arrondie la matrice_sucre avec 4 décimale
    new_mat = np.round(new_mat,4)
    return new_mat




def generate_gaussian_matrix(size, fwhm=3, amplitude=1):
    """
    Crée une matrice_sucre 2D avec une distribution gaussienne centrée.
    
    - size : taille de la matrice_sucre (carrée)
    - fwhm : largeur à mi-hauteur (Full Width Half Maximum)
    - amplitude : hauteur du pic
    """
    x = np.arange(0, size, 1)
    y = x[:, np.newaxis]
    
    x0 = y0 = size // 2  # centre de la gaussienne
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))  # conversion FWHM → sigma
    
    gaussian = amplitude * np.exp(-((x - x0)**2 + (y - y0)**2) / (2 * sigma**2))
    return gaussian
"""
# Paramètres
size = 81
fwhm = 30
amplitude = 1

# Génération et affichage
mat = generate_gaussian_matrix(size, fwhm, amplitude)

plt.imshow(mat, cmap='viridis', origin='lower')
plt.colorbar(label="Intensité gaussienne")
plt.title("Distribution gaussienne centrée dans une matrice_sucre")
plt.show()

"""