from fonctionNormeAngleToComp import compVect
from valeurs import constantes as cte
from Trajectoire import derivee_solve_ivp, distance
from Calcul_trajectoire import euler
from RechercheRacine import secante, bissection, recursive

from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

import time
from threading import Thread

constantes = cte # Permet de modifier les valeurs dans angle

def minimumDistanceVaisseau(angleAst, angleVaisseau):
    # print("Début de minimumDistanceVaisseau")
    global constantes

    # start = time.time()
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
        cond_init[2][i] = compVect(angleVaisseau, 6000/1e9)[i] + constantes['valeursInit'][0][i]

    # print("Init : " + str(time.time() - start))

#    traj = euler(constantes['pas'], constantes['intervalleSimulation'], cond_init)
    # start = time.time()
    sol = solve_ivp(derivee_solve_ivp, constantes['intervalleSimulation'], np.reshape(cond_init, nbAstres*5), rtol=constantes['r-tol'], atol=constantes['a-tol'], max_step=constantes['pas'])
    # print("Solve ivp : " + str(time.time() - start))

    # Traduction de sol dans traj
    # start = time.time()
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
    
    # print("End : " + str(time.time() - start))

    constantes = cte

    return float(answer)

def angleVaisseau(posInitTerre, vitInitTerre, posInitAst, normeVitInitAst, angleAst, normeVitInitVaisseau):
    global constantes

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
    angles = np.arange(0, 2*np.pi, (2*np.pi)/constantes['nbAngles'])
    answer = np.zeros_like(angles)

    start = time.time()
    for i in range(len(angles)):
        answer[i] = minimumDistanceVaisseau(angleAst, angles[i])
    print(answer)
    print("Time bornes : " + str(time.time() - start))
    
    start = time.time()
    reponse, state = recursive(lambda angleVaisseau : minimumDistanceVaisseau(angleAst, angleVaisseau), 0, 2*np.pi, constantes['tol'], constantes['tol_bornes'], 0.5, constantes['nbAngles'])
    print("Temps de recursive : " + str(time.time() - start))
    print("State : " + str(state))
    
    constantes = cte

    return reponse

def collision(angleVaisseau, angleAst):
    # print("Début de minimumDistanceVaisseau")
    global constantes

    # start = time.time()
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
        cond_init[2][i] = compVect(angleVaisseau, 6000/1e9)[i] + constantes['valeursInit'][0][i]

    # print("Init : " + str(time.time() - start))

#    traj = euler(constantes['pas'], constantes['intervalleSimulation'], cond_init)
    # start = time.time()
    sol = solve_ivp(derivee_solve_ivp, constantes['intervalleSimulation'], np.reshape(cond_init, nbAstres*5), rtol=constantes['r-tol'], atol=constantes['a-tol'], max_step=constantes['pas'])
    # print("Solve ivp : " + str(time.time() - start))

    # Traduction de sol dans traj
    # start = time.time()
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
    
    found = False
    i = 0
    while not found and i < len(points_x[0]):
        if distance(points_x[1][i], points_y[1][i], points_x[2][i], points_y[2][i]) < constantes['tol']:
            found = True
        i += 1
    if found == False:
        print("L'astéroïde n'entre pas en collision avec le vaisseau")
        return -1

    pas = constantes['intervalleSimulation'][1]/len(points_x[0])
    intervalleNouveau = [0, i*pas]
    nouvelleCondInit = traj[i]
    nouvelleCondInit[1][2] = (nouvelleCondInit[1][4]*nouvelleCondInit[1][2] + nouvelleCondInit[2][4]*nouvelleCondInit[2][2])/(nouvelleCondInit[1][4] + nouvelleCondInit[2][4])
    nouvelleCondInit[1][4] += nouvelleCondInit[2][4]
    nouvelleCondInit = nouvelleCondInit[:-1]

    sol = solve_ivp(derivee_solve_ivp, intervalleNouveau, np.reshape(nouvelleCondInit, nbAstres*5), rtol=constantes['r-tol'], atol=constantes['a-tol'], max_step=constantes['pas'])

    # Traduction de sol dans traj
    traj = np.zeros((len(sol.t), nbAstres, 5))
    for n in range(len(traj)):
        for j in range(nbAstres):
            for k in range(5):
                traj[n][j][k] = sol.y[j*5 + k][n]
    
    points_x_new = np.zeros((nbAstres, len(traj)))
    points_y_new = np.zeros((nbAstres, len(traj)))
    for j in range(len(traj)):
        for n in range(nbAstres):
            points_x_new[n][j] = traj[j][n][0]
            points_y_new[n][j] = traj[j][n][1]
    



# Fonction qui sert à signer l'angle avec la dérivée
def minimumDistancePourAngleVaisseau(angleVaisseau, angleAst, derivee, angles):
    gardien = True
    indice = -1
    for i in range(len(derivee)):
        if gardien and angleVaisseau >= angles[i] and angleVaisseau < angles[i+1]:
            indice = i
            gardien = False
    
    return minimumDistanceVaisseau(angleAst, angleVaisseau) * derivee[indice]/np.abs(derivee[indice])