import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random as rd
from repartition_sucre import *
from bacteria_mvt import GLUCOSE_CONSUMPTION,Bacteria
#from Energie_bact import *

# ===================== PARAMÈTRES =====================
M_TAILLE = 20                     # Taille de la matrice
NB_BACTERIES_INIT = 1       # Nombre de bactéries au début
ADD_GLUCOSE_INTERVAL = 10       # (non utilisé ici, mais prévu)
DUREE_GLUCOSE = 20        # Jusqu'à quand on ajoute du sucre
DUREE_DIFFUSION = 10000        # Jusqu'à quand on diffuse le sucre
DEBUT_BACTERIES = 1            # Quand les bactéries commencent à agir
MAX_ITER = 1000                # Pour un test rapide (change selon besoin)
SUCRE_POS = (M_TAILLE//2, M_TAILLE//2)              # Position d'injection de sucre
SUCRE_CONCENTRATION = 0.04      # Concentration injectée
SUCRE_RAYON = M_TAILLE//2                 # Rayon de diffusion du sucre
BACTERIE_COLOR = 'white'
BACTERIE_SIZE = 10
VITESSE_DIFFUSTION_SUCRE = 3  # Vitesse de diffusion du sucre

# ===================== INITIALISATION =====================
mat = init_matrice(M_TAILLE)
list_b = [Bacteria(rd.uniform(0, 1), rd.uniform(0, 1)) for _ in range(NB_BACTERIES_INIT)]

# Setup animation
fig, ax = plt.subplots()
im = ax.imshow(mat, cmap='turbo', vmin=0, vmax=1, extent=[0, M_TAILLE, 0, M_TAILLE], origin='lower')
plt.colorbar(im, ax=ax)
title = ax.set_title("Bioreacteur avec sucre et bactéries")
x_bact = [b.posx * M_TAILLE for b in list_b]
y_bact = [b.posy * M_TAILLE for b in list_b]
sc = ax.scatter(x_bact, y_bact, c=BACTERIE_COLOR, s=BACTERIE_SIZE)

# ===================== FONCTION D'ANIMATION =====================
def update(frame):
    global mat, list_b
    steps_per_frame = 1
    for _ in range(steps_per_frame):
        iteration = frame * steps_per_frame + _

        # --- Ajout de glucose ---
        if iteration < DUREE_GLUCOSE:
            mat = sucre_input(mat, SUCRE_POS, SUCRE_CONCENTRATION, SUCRE_RAYON)

        # --- Diffusion (loi de Fick) ---
        if iteration < DUREE_DIFFUSION:
            for i in range(VITESSE_DIFFUSTION_SUCRE):
                mat = update_sucre(mat)

        # --- Bactéries (après un certain temps) ---
        if iteration > DEBUT_BACTERIES:
            new_bacts = []
            for bacterie in list_b[:]:  # Copie pour éviter modification pendant itération
                if bacterie.death:
                    list_b.remove(bacterie)
                    continue
                bacterie.update_b_pos(mat)
                bacterie.update_state(mat)
                new_bact = bacterie.update_death_and_mitosis(mat, list_b)
                if isinstance(new_bact, Bacteria):
                    new_bacts.append(new_bact)
                mat = bacterie.update_eat(mat,GLUCOSE_CONSUMPTION)
            list_b.extend(new_bacts)

    # --- Affichage / mise à jour graphique ---
    x_bact = [b.posx * M_TAILLE for b in list_b]
    y_bact = [b.posy * M_TAILLE for b in list_b]
    colors = ['violet' if b.consommation_state == 'fermentation' else 'white' for b in list_b]
    sc.set_offsets(np.c_[x_bact, y_bact])
    sc.set_color(colors)
    im.set_array(mat)

    # --- Affichage console ---
    print(len(list_b), "bactéries vivantes")
    print("Nb bact fermentation :", len([b for b in list_b if b.consommation_state == 'fermentation']))
    print("Nb bact respiration :", len([b for b in list_b if b.consommation_state == 'respiration']))
    if iteration == MAX_ITER:
        print("Fin de la simulation.")
        exit()

    return [im, sc]

# ===================== LANCEMENT =====================
ani = animation.FuncAnimation(fig, update, frames=10_000, interval=100, blit=True)
plt.show()
