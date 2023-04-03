import numpy as np
import matplotlib.pyplot as plt
import json

from valeurs import constantes as cte

from fonctionNormeAngleToComp import compVect
from scipy.integrate import solve_ivp

import time

constantes = cte

def deriver(t, positions_vitesses_masse):
    nbAstres = len(positions_vitesses_masse)

    i = 0
    while i < nbAstres:
        to_check = [k for k in range(i+1, nbAstres)]
        for j in to_check:
            if positions_vitesses_masse[i][0] == positions_vitesses_masse[j][0] and positions_vitesses_masse[i][1] == positions_vitesses_masse[j][1]:
                return -1
        i += 1

    # Création du tableau renvoyé
    derivees = np.zeros((nbAstres, 5))

    for i in range(nbAstres):
        derivees[i][0] = positions_vitesses_masse[i][2]
        derivees[i][1] = positions_vitesses_masse[i][3]

    for i in range(nbAstres):
        ax = 0
        ay = 0

        # x cos, y sin
        # Accélération due au soleil
        if (positions_vitesses_masse[i][0]**2 + positions_vitesses_masse[i][1]**2 == 0):
            return -1
        ax += (constantes['cavendish'] * constantes['masseSoleil'] * (positions_vitesses_masse[i][0] / distance(positions_vitesses_masse[i][0], positions_vitesses_masse[i][1], 0, 0)))/ \
            (positions_vitesses_masse[i][0]**2 + positions_vitesses_masse[i][1]**2)
        ay += (constantes['cavendish'] * constantes['masseSoleil'] * (positions_vitesses_masse[i][1] / distance(positions_vitesses_masse[i][0], positions_vitesses_masse[i][1], 0, 0)))/ \
            (positions_vitesses_masse[i][0]**2 + positions_vitesses_masse[i][1]**2)
    
        # Calcul de l'accélération de l'astre par rapport à tous les autres du tableau
        autresAstres = [k for k in range(nbAstres)]
        autresAstres.remove(i)
        for j in autresAstres:
            # print(distance(positions_vitesses_masse[i][0], positions_vitesses_masse[i][1], positions_vitesses_masse[j][0], positions_vitesses_masse[j][1]))
            ax += (constantes['cavendish'] * positions_vitesses_masse[j][4] * ((positions_vitesses_masse[i][0] - positions_vitesses_masse[j][0]) / distance(positions_vitesses_masse[i][0], positions_vitesses_masse[i][1], positions_vitesses_masse[j][0], positions_vitesses_masse[j][1])))/ \
                (distance(positions_vitesses_masse[i][0], positions_vitesses_masse[i][1], positions_vitesses_masse[j][0], positions_vitesses_masse[j][1])**2)
            ay += (constantes['cavendish'] * positions_vitesses_masse[j][4] * ((positions_vitesses_masse[i][1] - positions_vitesses_masse[j][1]) / distance(positions_vitesses_masse[i][0], positions_vitesses_masse[i][1], positions_vitesses_masse[j][0], positions_vitesses_masse[j][1])))/ \
                (distance(positions_vitesses_masse[i][0], positions_vitesses_masse[i][1], positions_vitesses_masse[j][0], positions_vitesses_masse[j][1])**2)
        
        # Filtre les valeurs aberrantes pour le vaisseau quand il est trop proche de la Terre
        if i == 2 and (ax**2 + ay**2) > (1e8)**2:
            ax = 0
            ay = 0

        # L'accélération est vers les autres astres
        derivees[i][2] = -1 * ax
        derivees[i][3] = -1 * ay
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
def simulation(angleAst, angleVaisseau):
    if len(constantes['valeursInit']) < 3:
        constantes['valeursInit'].append([0,0,0,0, constantes['masseVaisseau']]) # Rajout du vaisseau
    nbAstres = len(constantes['valeursInit'])

    # Vitesse de l'astéroide
    vx, vy = compVect(angleAst, constantes['vitesseAsteroide'])
    cond_init = constantes['valeursInit']
    cond_init[1][2] = vx
    cond_init[1][3] = vy

    # Vitesse du vaisseau
    vx, vy = compVect(angleVaisseau, constantes['vitesseVaisseau'])
    cond_init[2][2] = vx + cond_init[0][2]
    cond_init[2][3] = vy + cond_init[0][3]

    # Position du vaisseau
    for i in [0,1]:
        cond_init[2][i] = compVect(angleVaisseau, 6371/1e9)[i] + constantes['valeursInit'][0][i]

#    traj = euler(constantes['pas'], constantes['intervalleSimulation'], constantes['valeursInit'])
    # start = time.time()
    sol = solve_ivp(derivee_solve_ivp, constantes['intervalleSimulation'], np.reshape(cond_init, nbAstres*5), rtol=constantes['r-tol'], atol=constantes['r-tol'], max_step=constantes['pas'])
    # print("Time : " + str(time.time() - start))

    # Traduction de sol dans traj
    traj = np.zeros((len(sol.t), nbAstres, 5))
    for i in range(len(traj)):
        for j in range(nbAstres):
            for k in range(5):
                traj[i][j][k] = sol.y[j*5 + k][i]

    # Extracting points coordinates
    nbAstres = len(traj[0])
    points_x = np.zeros((nbAstres, len(traj)))
    points_y = np.zeros((nbAstres, len(traj)))
    # print("Making points arrays")
    for i in range(len(traj)):
        for j in range(nbAstres):
            points_x[j][i] = traj[i][j][0]
            points_y[j][i] = traj[i][j][1]

    # Exporting to file
    # json_points = json.dumps([points_x, points_y], indent=4)
    # with open("points_" + str(constantes['pas']) + "_" + str(constantes['intervalleSimulation']) + "_.json", "w") as outfile:
    #     outfile.write(json_points)

    # Plotting the trajectory
    plt.axis('equal')

    plt.plot(points_x[0], points_y[0], 'green')
    plt.plot(points_x[1], points_y[1], 'brown')
    plt.plot(points_x[2], points_y[2], 'grey')
    plt.scatter(0, 0, color='orange')
    plt.scatter(cte['valeursInit'][0][0], cte['valeursInit'][0][1], color='black')

    plt.title("Space simulation")
    plt.xlabel("X coordinate")
    plt.ylabel("Y coordinate")
    plt.legend(["Earth", "Asteroid", "Spaceship", "Sun"], loc='best')

    plt.show()
    
def derivee_solve_ivp(t, positions_vitesses_masses):
    nbAstre = int(len(positions_vitesses_masses)/5)

    arr = np.zeros((nbAstre, 5))
    for i in range(nbAstre):
        for j in range(5):
            arr[i][j] = positions_vitesses_masses[i*5+j]

    new_arr = deriver(t, arr)
    if type(new_arr) is int and new_arr == -1:
        new_arr = arr
        for i in range(nbAstre):
            new_arr[i][0] += t*new_arr[i][2]
            new_arr[i][1] += t*new_arr[i][3]
    new_arr = np.reshape(new_arr, nbAstre*5)
    return new_arr

def distance(x1, y1, x2, y2):
    return np.sqrt((x1 - x2)**2 + (y1 - y2)**2)