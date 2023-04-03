from fonctionNormeAngleToComp import compVect
from valeurs import constantes as cte
from Calcul_trajectoire import euler, deriver, derivee_solve_ivp, distance
from RechercheRacine import secante

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

constantes = cte # Permet de modifier les valeurs dans angle

def minimumDistance(angle):
    # print("Début de minimumDistance")

    nbAstres = len(constantes['valeursInit'])

    vx, vy = compVect(angle, constantes['vitesseAsteroide'])
    cond_init = constantes['valeursInit']
    cond_init[1][2] = vx
    cond_init[1][3] = vy

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

    answer = distance(points_x[0][0], points_y[0][0], points_x[1][0], points_y[1][0])
    for i in range(len(traj)):
        distancePotentielle = distance(points_x[0][i], points_y[0][i], points_x[1][i], points_y[1][i])
        if distancePotentielle < answer:
            answer = distancePotentielle

    return float(answer)

def angle(posInitTerre, vitInitTerre, posInitAst, normeVitInitAst):
    # Insertion des paramètres dans constantes
    constantes['valeursInit'][0][0] = posInitTerre[0]
    constantes['valeursInit'][0][1] = posInitTerre[1]
    constantes['valeursInit'][0][2] = vitInitTerre[0]
    constantes['valeursInit'][0][3] = vitInitTerre[1]
    constantes['valeursInit'][1][0] = posInitAst[0]
    constantes['valeursInit'][1][1] = posInitAst[1]
    constantes['vitesseAsteroide'] = normeVitInitAst

    # Création des tableaux pour la dérivée et les bornes 
    angles = np.arange(0, 2*np.pi, (2*np.pi)/20)
    answer = np.zeros_like(angles)

    for i in range(len(angles)):
        answer[i] = minimumDistance(angles[i])

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

    reponse, state = secante(lambda angle: minimumDistancePourAngle(angle, derivee, angles), borne_1, borne_2, constantes['tol'])

    return reponse

# Fonction qui sert à signer l'angle avec la dérivée
def minimumDistancePourAngle(angle, derivee, angles):
    gardien = True
    indice = -1
    for i in range(len(derivee)):
        if gardien and angle >= angles[i] and angle < angles[i+1]:
            indice = i
            gardien = False
    
    return minimumDistance(angle) * derivee[indice]/np.abs(derivee[indice])