import numpy as np
import matplotlib.pyplot as plt
import random
import pygame as pg


m_taille = 3 # Taille de la matrice

def init_matrice(taille):
    """Initialise une matrice de taille donnée avec des zeroes."""
    return np.zeros((taille, taille))


def update_sucre(old_mat):
    diffusion_mat=np.zeros((m_taille, m_taille))
    for i in range(old_mat.shape[0]-1):
        for j in range(old_mat.shape[1]-1):
            nb_de_voisin = 4
            for a in (1,-1):
                diffusion=old_mat[i][j]*(1-old_mat[i+a][j])/nb_de_voisin
                diffusion_mat[i+a][j] += diffusion
                print("avant",old_mat[i][j])
                diffusion_mat[i][j] -= diffusion
                print("diffusion",diffusion)
    print("diffusion",diffusion_mat)
    new_mat=old_mat +diffusion_mat
    return new_mat
            
            
mat = init_matrice(m_taille)
mat[1][1]+=1
print(mat)
print(mat.shape)
for i in range(3):
    mat=update_sucre(mat)
    print("new_mat",mat)