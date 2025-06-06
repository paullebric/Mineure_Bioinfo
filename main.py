import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random as rd
from repartition_sucre import *
from bacteria_mvt import *
from graph import *
from toxine import *
#from Energie_bact import *

# ===================== PARAMÈTRES =====================
SUCRE_BASE = 0 # Concentration de sucre de base dans la matrice_sucre
M_TAILLE = 25                    # Taille de la matrice_sucre
NB_BACTERIES_INIT = 3      # Nombre de bactéries au début
DUREE_GLUCOSE = 1000000        # Jusqu'à quand on ajoute du sucre
DUREE_DIFFUSION = 150000000       # Jusqu'à quand on diffuse le sucre
DEBUT_BACTERIES = 1            # Quand les bactéries commencent à agir
MAX_ITER = 20000           # Pour un test rapide (change selon besoin)
SUCRE_POS = (M_TAILLE//2, M_TAILLE//2)              # Position d'injection de sucre
SUCRE_CONCENTRATION = 0.2    # Concentration injectée
SUCRE_RAYON = 1                # Rayon de diffusion du sucre
BACTERIE_COLOR = 'white'
BACTERIE_SIZE = 5
VITESSE_DIFFUSTION_SUCRE = 3  # Vitesse de diffusion du sucre

DEATH = True  # Si True, les bactéries peuvent mourir
MITOSE = True  # Si True, les bactéries peuvent se diviser
EAT = True  # Si True, les bactéries peuvent manger
STATE = True  # Si True, les bactéries changent d'état de consommation
MOOVE = True
CONCERTE = False  # Si True, les bactéries se déplacent en concerté
# ===================== INITIALISATION =====================
mat = init_matrice_sucre(M_TAILLE, SUCRE_BASE)
mat_tox = init_mat_tox(M_TAILLE)

list_b = [Bacteria(rd.uniform(0, 1), rd.uniform(0, 1)) for _ in range(NB_BACTERIES_INIT)]

# Setup animation
fig, ax = plt.subplots()
im = ax.imshow(mat, cmap='turbo', vmin=0, vmax=1, extent=[0, M_TAILLE, 0, M_TAILLE], origin='lower')
plt.colorbar(im, ax=ax)
title = ax.set_title("Bioreacteur avec sucre et bactéries")
x_bact = [b.posx * M_TAILLE for b in list_b]
y_bact = [b.posy * M_TAILLE for b in list_b]
sc = ax.scatter(x_bact, y_bact, c=BACTERIE_COLOR, s=BACTERIE_SIZE)

bact_counts = []  # Pour stocker le nombre de bactéries à chaque itération
bact_counts_respi = []  # Pour les bactéries en respiration
bact_counts_ferment = []  # Pour les bactéries en fermentation
gluc_sum = []
tox_sum = []
# ===================== FONCTION D'ANIMATION =====================
def update(frame):
    global mat, list_b, mat_tox
    steps_per_frame = 10
    for _ in range(steps_per_frame):
        iteration = frame * steps_per_frame + _

        # --- Ajout de glucose ---
        if iteration < DUREE_GLUCOSE:
            mat = sucre_input(mat, SUCRE_POS, SUCRE_CONCENTRATION, SUCRE_RAYON)
            #mat = sucre_input(mat, (4,3), SUCRE_CONCENTRATION*1.2, SUCRE_RAYON)

        
        # --- Diffusion (loi de Fick) ---
        if iteration < DUREE_DIFFUSION:
            for i in range(VITESSE_DIFFUSTION_SUCRE):
                mat = update_sucre(mat)

        
        # --- Bactéries (après un certain temps) ---
        if iteration > DEBUT_BACTERIES:
            new_bacts = []
            for bacterie in list_b[:]:  # Copie pour éviter modification pendant itération
                
                if DEATH :
                    if bacterie.death:
                        list_b.remove(bacterie)
                        continue

                if MOOVE : 
                    if CONCERTE:
                        bacterie.update_Rij3(mat)
                    bacterie.update_b_pos(mat,CONCERTE,list_b)
                
                if STATE:bacterie.update_state(mat)

                if EAT:

                    matdouble = bacterie.update_eat(mat,mat_tox,GLUCOSE_CONSUMPTION)
                    mat = matdouble[0]
                    mat_tox = matdouble[1]
                
                new_bact = bacterie.update_death_and_mitosis(mat_tox, list_b)
                if MITOSE :
                    if isinstance(new_bact, Bacteria):
                        new_bacts.append(new_bact)
            list_b.extend(new_bacts)
            bact_counts.append(len(list_b))
            bact_counts_respi.append(len([b for b in list_b if b.consommation_state == 'respiration']))
            bact_counts_ferment.append(len([b for b in list_b if b.consommation_state == 'fermentation']))
            gluc_sum.append(np.sum(mat))
            tox_sum.append(np.sum(mat_tox))

    # --- Affichage / mise à jour graphique ---
    x_bact = [b.posx * M_TAILLE for b in list_b]
    y_bact = [b.posy * M_TAILLE for b in list_b]
    colors = ['violet' if b.consommation_state == 'fermentation' else 'white' for b in list_b]
    sc.set_offsets(np.c_[x_bact, y_bact])
    sc.set_color(colors)
    im.set_array(mat)

    # --- Affichage console ---
    
    


    if iteration > MAX_ITER or len(list_b)==0:

        ani.event_source.stop()
        afficher_croissance_bact(bact_counts)
        afficher_croissance_bact_by_state(bact_counts_respi,bact_counts_ferment,gluc_sum,tox_sum)

    return [im, sc]

# ===================== LANCEMENT =====================
ani = animation.FuncAnimation(fig, update, frames=10_000, interval=100, blit=True)
plt.show()
