import numpy as np
import matplotlib.pyplot as plt
import random as rd
from repartition_sucre import *
from bacteria_mvt import Bacteria

m_taille = 50 # Taille de la matrice
add_glucose = 10 # intervalle d'ajout de glucose

# Initialisation de la matrice et des bactéries
mat = init_matrice(m_taille)
list_b = [Bacteria(rd.uniform(0, 1), rd.uniform(0, 1)) for _ in range(1000)]

# Animation setup
fig, ax = plt.subplots()
im = ax.imshow(mat, cmap='turbo', vmin=0, vmax=1, extent=[0, m_taille, 0, m_taille], origin='lower')
plt.colorbar(im, ax=ax)
title = ax.set_title("Bioreacteur avec sucre et bactéries")

# Initialisation du scatter plot pour les bactéries
x_bact = [b.posx * m_taille for b in list_b]
y_bact = [b.posy * m_taille for b in list_b]
sc = ax.scatter(x_bact, y_bact, c='white', s=10)  # points blancs pour les bactéries

# Fonction d'animation = boucle principale
def update(frame):
    global mat
    steps_per_frame = 100  # vitesse de l'animation
    for _ in range(steps_per_frame):
        iteration = frame * steps_per_frame + _
        if iteration % add_glucose == 0:  # Ajout de glucose à chaque intervalle n défini
            if iteration < 200: # Ajout de glucose au début de la simulation jusqqquaa n simulaation
                #mat = sucre_input(mat, (m_taille // 2, m_taille // 2), 0.1, 1)
                mat = sucre_input(mat, (12, 10), 0.1, 1)
                mat = sucre_input(mat, (20, 22), 0.5, 1)
            if iteration < 400: # Loi de flick = repartition du glucose jusqu'à n frames
                mat = update_sucre(mat)
        #les bacteries existe a partir de ce moment
        if iteration> 1:
            for bacterie in list_b:
                if bacterie.death == True:
                    list_b.remove(bacterie)
                # Mise à jour de la position des bactéries
                bacterie.update_b_pos(mat)
                # Mise à jour de la mort et de la mitose des bactéries
                if False :
                    new_bact = bacterie.update_death_and_mitosis(mat,list_b)
                    if isinstance(new_bact,Bacteria):
                        list_b.append(new_bact)
                    mat = bacterie.update_eat(mat)
    im.set_array(mat)
    # Mise à jour des positions des bactéries
    x_bact = [b.posx * m_taille for b in list_b]
    y_bact = [b.posy * m_taille for b in list_b]
    sc.set_offsets(np.c_[x_bact, y_bact])
    return [im, sc]

ani = animation.FuncAnimation(fig, update, frames=10, interval=100, blit=True)
plt.show()



