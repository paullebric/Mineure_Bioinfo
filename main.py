import numpy as np
import matplotlib.pyplot as plt
import random
from repartition_sucre import *
from bacteria_mvt import *

mat = init_matrice(m_taille)

list_b = [Bacteria(random.uniform(0, 1), random.uniform(0, 1)) for _ in range(100)]

# Animation setup
fig, ax = plt.subplots()
im = ax.imshow(mat, cmap='turbo', vmin=0, vmax=1)
plt.colorbar(im, ax=ax)
title = ax.set_title("Diffusion de sucre")

# Initialisation du scatter plot pour les bactéries
x_bact = [b.posx * m_taille for b in list_b]
y_bact = [b.posy * m_taille for b in list_b]
sc = ax.scatter(x_bact, y_bact, c='white', s=10)  # points blancs pour les bactéries

# Fonction d'animation
def update(frame):
    global mat
    mat = update_sucre(mat)
    if frame % add_glucose == 0:
        mat = sucre_input(mat, (m_taille // 2, m_taille // 2), 0.5, 1)
    for bacterie in list_b:
        bacterie.update_b_pos(mat)
    im.set_array(mat)
    # Mise à jour des positions des bactéries
    x_bact = [b.posx * m_taille for b in list_b]
    y_bact = [b.posy * m_taille for b in list_b]
    sc.set_offsets(np.c_[x_bact, y_bact])
    return [im, sc]

ani = animation.FuncAnimation(fig, update, frames=10000, interval=0.1, blit=True)
plt.show()



