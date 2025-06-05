import numpy as np
import matplotlib.pyplot as plt
import random as rd
from repartition_sucre import *
from bacteria_mvt import Bacteria
from Energie_bact import *
m_taille = 3 # Taille de la matrice
add_glucose = 10 # intervalle d'ajout de glucose

# Initialisation de la matrice et des bactéries
mat = init_matrice(m_taille)
list_b = [Bacteria(rd.uniform(0, 1), rd.uniform(0, 1)) for _ in range(3)]

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
    steps_per_frame = 1  # vitesse de l'animation
    for _ in range(steps_per_frame):
        #gradglobal(list_b, mat)  # Mise à jour du gradient global pour chaque bactérie
        iteration = frame * steps_per_frame + _
        if iteration < 200000: # Ajout de glucose au début de la simulation jusqqquaa n simulaation
                mat = sucre_input(mat, (0, 0), 0.3, 1)
                # mat = sucre_input(mat, (12, 10), 0.5, 1)
                # mat = sucre_input(mat, (20, 22), 1, 1)
        if iteration < 10000: # Loi de flick = repartition du glucose jusqu'à n frames
                mat = update_sucre(mat)
        #les bacteries existe a partir de ce moment
        if True:
            if iteration > 20:
                for bacterie in list_b:
                    if bacterie.death == True:
                        list_b.remove(bacterie)
                    # Mise à jour de la position des bactéries
                    bacterie.update_b_pos(mat)
                    if True :
                        new_bact = bacterie.update_death_and_mitosis(mat,list_b)
                        if isinstance(new_bact,Bacteria):
                            list_b.append(new_bact)
                        mat = bacterie.update_eat(mat)
    print(len(list_b), "bactéries vivantes")
    for bact in list_b:
        print(bact.gradnumx, bact.gradnumy)
    if iteration == 1000:print("fin"),exit()
    im.set_array(mat)
    # Mise à jour des positions des bactéries
    x_bact = [b.posx * m_taille for b in list_b]
    y_bact = [b.posy * m_taille for b in list_b]
    sc.set_offsets(np.c_[x_bact, y_bact])
    return [im, sc]

ani = animation.FuncAnimation(fig, update, frames=10000, interval=100, blit=True)
plt.show()



