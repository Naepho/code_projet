import numpy as np
import matplotlib.pyplot as plt
import json

from valeurs import constantes

def deriver(t, positions_vitesses_masse):
    nbAstres = len(positions_vitesses_masse)

    i = 0
    while i < nbAstres:
        to_check = [k for k in range(i+1, nbAstres)]
        for j in to_check:
            if positions_vitesses_masse[i][0] == positions_vitesses_masse[j][0] and positions_vitesses_masse[i][1] == positions_vitesses_masse[j][1]:
                return -1
        i += 1

    # création du tableau renvoyé
    derivees = np.zeros((nbAstres, 5))

    for i in range(nbAstres):
        derivees[i][0] = positions_vitesses_masse[i][2]
        derivees[i][1] = positions_vitesses_masse[i][3]

    for i in range(nbAstres):
        ax = 0
        ay = 0

        # x cos, y sin
        # Accélération due au soleil
        ax -= (constantes['cavendish'] * constantes['masseSoleil'] * (positions_vitesses_masse[i][0] / np.sqrt(positions_vitesses_masse[i][0]**2 + positions_vitesses_masse[i][1]**2)))/ \
            (positions_vitesses_masse[i][0]**2 + positions_vitesses_masse[i][1]**2)
        ay -= (constantes['cavendish'] * constantes['masseSoleil'] * (positions_vitesses_masse[i][1] / np.sqrt(positions_vitesses_masse[i][0]**2 + positions_vitesses_masse[i][1]**2)))/ \
            (positions_vitesses_masse[i][0]**2 + positions_vitesses_masse[i][1]**2)
    
        # Calcul de l'accélération de l'astre par rapport à tous les autres du tableau
        autresAstres = [k for k in range(nbAstres)]
        autresAstres.remove(i)
        for j in autresAstres:
            ax -= (constantes['cavendish'] * positions_vitesses_masse[j][4] * ((positions_vitesses_masse[i][0] - positions_vitesses_masse[j][0]) / np.sqrt((positions_vitesses_masse[i][0] - positions_vitesses_masse[j][0])**2 + (positions_vitesses_masse[i][1] - positions_vitesses_masse[j][1])**2)))/ \
                ((positions_vitesses_masse[i][0] - positions_vitesses_masse[j][0])**2 + (positions_vitesses_masse[i][1] - positions_vitesses_masse[j][1])**2)
            ay -= (constantes['cavendish'] * positions_vitesses_masse[j][4] * ((positions_vitesses_masse[i][1] - positions_vitesses_masse[j][1]) / np.sqrt((positions_vitesses_masse[i][0] - positions_vitesses_masse[j][0])**2 + (positions_vitesses_masse[i][1] - positions_vitesses_masse[j][1])**2)))/ \
                ((positions_vitesses_masse[i][0] - positions_vitesses_masse[j][0])**2 + (positions_vitesses_masse[i][1] - positions_vitesses_masse[j][1])**2)

        # L'accélération est vers les autres astres
        derivees[i][2] = ax
        derivees[i][3] = ay
        derivees[i][4] = 0

    return derivees

def euler(pas, intervalle, cond_init):
    # print("Starting the Euler method")

    # Time array
    t = np.arange(intervalle[0], intervalle[1] + pas, pas)
    
    # Answer  array
    answer = np.zeros((len(t), len(cond_init), len(cond_init[0])))

    # Setting inital conditions
    for i in range(len(cond_init)):
        answer[0][i] = cond_init[i]

    # The main loop
    for i in range(1, len(t)):
        answer[i] = answer[i-1] + pas*deriver(t[i], answer[i-1])
    
    return answer

# Calcule juste les trajectoires et les affiche
def simulation():
    traj = euler(constantes['pas'], constantes['intervalleSimulation'], constantes['valeursInit'])

    # Extracting points coordinates
    nbAstres = len(traj[0])
    points_x = np.zeros((nbAstres, len(traj)))
    points_y = np.zeros((nbAstres, len(traj)))
    print("Making points arrays")
    for i in range(len(traj)):
        for j in range(nbAstres):
            points_x[j][i] = traj[i][j][0]
            points_y[j][i] = traj[i][j][1]

    # Exporting to file
    json_points = json.dumps([points_x, points_y], indent=4)
    with open("points_" + str(constantes['pas']) + "_" + str(constantes['intervalleSimulation']) + "_.json", "w") as outfile:
        outfile.write(json_points)

    # Plotting the trajectory
    plt.axis('equal')

    plt.plot(points_x[0], points_y[0], 'g')
    plt.plot(points_x[1], points_y[1], 'brown')
    plt.scatter(0, 0, color='orange')

    plt.title("Space simulation")
    plt.xlabel("X coordinate")
    plt.ylabel("Y coordinate")
    plt.legend(["Earth", "Asteroid", "Sun"], loc='best')

    plt.show()