# Mineure_Bioinfo
Simulation of a Bioreactor

En utilisant numpy
une matrice pour la concentration en sucre + une matrice pour le nombre de bacteria (chaque case est une position dans le bioreacteur)

Déplacement = Csucre * 0.25

"""""" CODE """"""
****matrice concentration en sucre --> le sucre arrive par la case du milieu (constante d'arrivée)
	LE sucre se diffuse selon le gradient
	formule = k1(Cd-Ca)*((1-Cd)/4) avec Cd concentration de la case qu'on regarde et Ca concentration de la case d'arrivé
	Rajouter des entrées de glucose
	coefficients différents sur x et y pour faire une ellispe 
	2 types de grille (homogène et hétérogène)

	
****matrice des bactéries
	position de départ soit aléatoire (probabilité d'apparition pour chaque case) soit cluster
	--> on génère des bactéries au hasard
	déplacement à gérer plus tard
	Est-ce que on gère les bactéries une par une (si oui comment on gère leur position et comment on l'affiche) 
	---> on va faire bactéries = points
	Est-ce que on choisit une approche ou est-ce qu'on les combine 
	--> on fait selon chemotaxie pour l'instant
	--> Formule déplacement : x(t+1) = x(t)+(l*Fx)*r	
		x(t), la position
		l, une constante
		Fx, le gradient de sucre 
			Fx = pos_case*(Cd-Cg)
		r, un random 
	--> Même formule pour y
****Liste constantes : 
concentration sucre initiale
nombre de bactéries
death rate
mitosis rate
eat rate


****class des bactéries ==>
	nombre de mitose avant mort
	threshold de mitose
	qté de sucre bouffé par les bactéries ==> peut etre gradient de la capacité de bouffe de sucre par les bactéries.
	Nombres de vie (nombre de fois d'affilé ou la bactérie peut ne pas manger assez de sucre)
	bacteries.update : fonction qui update l'etat de la bacterie
	déplacement aléatoire impacté par la qtité de sucre
	
****A chaque iteration (temps qui avance) ==>
	une fonction qui gère :la repartition du sucre dans l'espace au fil du temps
	une fonction qui gère :la conso de sucre par les Bacteries
	une fonction qui gère :le deplacement des bactéries
	bacteries.update sur toute les bacteries
	
	une fonction qui affiche la matrice des bactéries avec en fond la concentration en sucre

	deplacement des bactéries ==> on normalise 2 fois
""""" FIN CODE """""




Project Main Components:
	Find a classic Bioreactor fonctionning to based our model on it (sugar repartition ect…) maybe the bioreactor fermentation
	Define the Space:
	2D grid

Define the Objects:
	Bacteria and Nutrients ==> choose a bacteria (build a graph of it's interactions with sugar/nutrients ect…)

Define Object Movement:
	bacteria mooves randomly with more chance to go in direction of the sugar gradient

Nutrients (Sugar):
	Follow Brownian motion or diffusion.
	Implement heterogeneous distribution initially.
	Maybe look at how the suger is distributed in different bioreactors

Bacteria:
	Growth, death, and state changes. ==> how they die the chance of them dying what affect their death chance ect… (gillepsie algorithm)
	Move according to Chemotaxis (follow sugar gradients).

Various movement models:
	Direct gradient-based motion.
	Monte Carlo + Simulated Annealing for global optimization of biomass production.

Visualization:
	2D graphical output showing bacteria and sugar distributions. (maybe with pygame)

