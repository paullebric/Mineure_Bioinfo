import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random
import pygame as pg


m_taille = 20 # Taille de la matrice
add_glucose = 50 # nb d'itérations avant l'ajout de glucose

def init_matrice(taille):
    """Initialise une matrice de taille donnée avec des zeroes."""
    return np.zeros((taille, taille))

def sucre_input(mat,position,amount,area):
    x=position[0]
    y=position[1]
    for i in range(x-area,x+area):
        for j in range(y-area,y+area):
            print("x",x,"y",y)
            if mat[i][j]+amount>1:
                mat[i][j]=1
            else:
                mat[i][j]+=amount
    return mat

def afficher_matrice(mat, titre=""):
    plt.imshow(mat, cmap='YlOrRd', vmin=0, vmax=1)
    plt.colorbar(label='Concentration de sucre')
    plt.title(titre)
    plt.axis('off')  # cache les axes
    plt.show()

def update_sucre(old_mat):
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    diffusion_mat=np.zeros((m_taille, m_taille))
    for i in range(m_taille):
        for j in range(m_taille):
            nb_de_voisin = 4
            for di, dj in directions:
                ni, nj = i + di, j + dj
                #update les 4 diffusions en sucre dans les cases voisines
                if 0 <= ni < m_taille and 0 <= nj < m_taille:
                    diffusion=(old_mat[i][j]-old_mat[ni][nj])*(1-old_mat[ni][nj])/nb_de_voisin
                    if diffusion>0:
                        diffusion_mat[ni][nj] += diffusion
                        diffusion_mat[i][j] -= diffusion
    new_mat=old_mat +diffusion_mat
    return new_mat
            
            
mat = init_matrice(m_taille)

#animation
fig, ax = plt.subplots()
im = ax.imshow(mat, cmap='YlOrRd', vmin=0, vmax=1)
plt.colorbar(im, ax=ax)
plt.title("Diffusion de sucre")
title = ax.set_title("Diffusion de sucre")
# Fonction d'animation
def update(frame):
    global mat
    mat = update_sucre(mat)
    if frame%add_glucose==0:
        mat = sucre_input(mat,(m_taille//2, m_taille//2), 0.5,1)
    im.set_array(mat)
    return [im]

# Lancement de l'animation
ani = animation.FuncAnimation(fig, update, frames=100000, interval=1, blit=True)
plt.show()
