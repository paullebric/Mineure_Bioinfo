import numpy as np
import matplotlib.pyplot as plt
import random
from repartition_sucre import *
from bacteria_mvt import *
from interaction_bg import *

m_taille = 5 # Taille de la matrice
add_glucose = 10 # intervalle d'ajout de glucose

# Initialisation de la matrice et des bactéries
mat = init_matrice(m_taille)
list_b = [Bacteria(random.uniform(0, 1), random.uniform(0, 1)) for _ in range(100)]

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
    if frame <4:
        mat = update_sucre(mat)
    if frame % add_glucose == 0:
        if frame <10:
            mat = sucre_input(mat, (m_taille // 2, m_taille // 2), 0.5, 1)
    for bacterie in list_b:
        if bacterie.death == True:
            list_b.remove(bacterie)
        bacterie.update_b_pos(mat)
        bacterie.update_death_and_mitosis(mat,list_b)
    im.set_array(mat)
    # Mise à jour des positions des bactéries
    x_bact = [b.posx * m_taille for b in list_b]
    y_bact = [b.posy * m_taille for b in list_b]
    sc.set_offsets(np.c_[x_bact, y_bact])
    return [im, sc]

ani = animation.FuncAnimation(fig, update, frames=10000, interval=1, blit=True)
plt.show()



