import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random as rd


def afficher_croissance_bact(bacteria_counts):
    plt.figure()
    plt.plot(bacteria_counts, label="Population bactérienne", color='green')
    plt.xlabel("Temps (frame)")
    plt.ylabel("Nombre de bactéries vivantes")
    plt.title("Croissance bactérienne au cours du temps")
    plt.legend()
    plt.grid(True)
    plt.show()

def afficher_croissance_bact_by_state(bacteria_counts_respi,bacteria_counts_ferment,gluc_sum,tox_sum):
    plt.figure()
    plt.plot(bacteria_counts_respi, label="Population bactérienne (respiration)", color='black')
    plt.plot(bacteria_counts_ferment, label="Population bactérienne (fermentation)", color='violet')
    plt.plot(gluc_sum, label="Concentration de glucose", color='orange')
    plt.plot(tox_sum, label="Concentration de toxine", color='red')
    plt.xlabel("Temps (frame)")
    plt.ylabel("Nombre de bactéries vivantes")
    plt.title("Croissance bactérienne au cours du temps")
    plt.legend()
    plt.grid(True)
    plt.show()