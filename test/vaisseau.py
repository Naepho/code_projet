from fonctionNormeAngleToComp import compVect
from valeurs import constantes as cte
from Trajectoire import derivee_solve_ivp, distance
from Calcul_trajectoire import euler
from RechercheRacine import secante

from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

import time

constantes = cte # Permet de modifier les valeurs dans angle

def minimumDistanceVaisseau(angleAst, angleVaisseau):
    # print("Début de minimumDistanceVaisseau")

    constantes['valeursInit'].append([0,0,0,0, constantes['masseVaisseau']]) # Rajout du vaisseau
    nbAstres = len(constantes['valeursInit'])

    # Vitesse de l'astéroide
    vx, vy = compVect(angleAst, constantes['vitesseAsteroide'])
    cond_init = constantes['valeursInit']
    cond_init[1][2] = vx
    cond_init[1][3] = vy

    # Vitesse du vaisseau
    vx, vy = compVect(angleVaisseau, constantes['vitesseVaisseau'])
    cond_init[2][2] = vx + cond_init[1][2]
    cond_init[2][3] = vy + cond_init[1][3]

    # Position du vaisseau
    for i in [0,1]:
        cond_init[2][i] = compVect(angleVaisseau, 6000/1e9)[i] + constantes['valeursInit'][0][i]

#    traj = euler(constantes['pas'], constantes['intervalleSimulation'], cond_init)
    sol = solve_ivp(derivee_solve_ivp, constantes['intervalleSimulation'], np.reshape(cond_init, nbAstres*5), rtol=constantes['r-tol'], atol=constantes['r-tol'], max_step=constantes['pas'])
    
    # Traduction de sol dans traj
    traj = np.zeros((len(sol.t), nbAstres, 5))
    for i in range(len(traj)):
        for j in range(nbAstres):
            for k in range(5):
                traj[i][j][k] = sol.y[j*5 + k][i]
    
    points_x = np.zeros((nbAstres, len(traj)))
    points_y = np.zeros((nbAstres, len(traj)))
    for i in range(len(traj)):
        for j in range(nbAstres):
            points_x[j][i] = traj[i][j][0]
            points_y[j][i] = traj[i][j][1]

    answer = distance(points_x[1][0], points_y[1][0], points_x[2][0], points_y[2][0])
    for i in range(len(traj)):
        distancePotentielle = distance(points_x[1][i], points_y[1][i], points_x[2][i], points_y[2][i])
        if distancePotentielle < answer:
            answer = distancePotentielle

    return float(answer)

def angleVaisseau(posInitTerre, vitInitTerre, posInitAst, normeVitInitAst, angleAst, normeVitInitVaisseau):
    # Insertion des paramètres dans constantes
    constantes['valeursInit'][0][0] = posInitTerre[0]
    constantes['valeursInit'][0][1] = posInitTerre[1]
    constantes['valeursInit'][0][2] = vitInitTerre[0]
    constantes['valeursInit'][0][3] = vitInitTerre[1]
    constantes['valeursInit'][1][0] = posInitAst[0]
    constantes['valeursInit'][1][1] = posInitAst[1]
    constantes['vitesseAsteroide'] = normeVitInitAst
    constantes['vitesseVaisseau'] = normeVitInitVaisseau

    # Création des tableaux pour la dérivée et les bornes 
    angles = np.arange(0, 2*np.pi, (2*np.pi)/20)
    answer = np.zeros_like(angles)

    for i in range(len(angles)):
        answer[i] = minimumDistanceVaisseau(angleAst, angles[i])
    print(answer)

    derivee = np.zeros((len(answer) - 1))
    for i in range(len(answer) - 1):
        derivee[i] = answer[i+1] - answer[i]
    
    # Calcul des bornes
    borne_1 = -1
    borne_2 = -1
    tol_bornes = constantes['tol_bornes']
    below_tol = False
    for i in range(len(answer)):
        if not below_tol:
            if answer[i] <= tol_bornes:
                borne_1 = angles[i]
                below_tol = True
        else:
            if answer[i] > tol_bornes:
                borne_2 = angles[i]
                below_tol = False

    if borne_1 == -1 and borne_2 == -1:
        return -1
    
    if borne_1 != -1 and borne_2 == -1:
        borne_2 = 2*np.pi

#    reponse, state = secante(lambda angleVaisseau: minimumDistancePourAngleVaisseau(angleVaisseau, angleAst, derivee, angles), borne_1, borne_2, constantes['tol'])
    reponse = -1
    state = -1
    print("State :" + str(state))
    print("stds")
    print(state)

    # Plotting the trajectory
    plt.axis('equal')

    plt.plot(angles, answer, 'g')

    plt.title("Space simulation")
    plt.xlabel("Angle")
    plt.ylabel("minimumDistanceVaisseau")

    plt.show()

    return reponse

# Fonction qui sert à signer l'angle avec la dérivée
def minimumDistancePourAngleVaisseau(angleVaisseau, angleAst, derivee, angles):
    gardien = True
    indice = -1
    for i in range(len(derivee)):
        if gardien and angleVaisseau >= angles[i] and angleVaisseau < angles[i+1]:
            indice = i
            gardien = False
    
    return minimumDistanceVaisseau(angleAst, angleVaisseau) * derivee[indice]/np.abs(derivee[indice])