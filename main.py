import numpy as np
import matplotlib.pyplot as plt
import random
from repartition_sucre import *
from bacteria_mvt import *

mat = init_matrice(m_taille)
#animation
fig, ax = plt.subplots()
im = ax.imshow(mat, cmap='turbo', vmin=0, vmax=0.5)
plt.colorbar(im, ax=ax)
plt.title("Diffusion de sucre")
title = ax.set_title("Diffusion de sucre")
# Fonction d'animation
def update(frame):
    global mat
    mat = update_sucre(mat)
    if frame%add_glucose==0:
        mat = sucre_input(mat,(m_taille//2, m_taille//2), 0.5,1)
    print("Frame:", frame)
    im.set_array(mat)
    return [im]
# Lancement de l'animation
ani = animation.FuncAnimation(fig, update, frames=100000, interval=1, blit=True)
plt.show()




