from fonctionNormeAngleToComp import compVect
from valeurs import constantes as cte
from Trajectoire import derivee_solve_ivp, distance
from Calcul_trajectoire import euler
from RechercheRacine import secante, bissection

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

    derivee = np.zeros((len(answer) - 1))
    for i in range(len(answer) - 1):
        derivee[i] = answer[i+1] - answer[i]
    
    # Calcul des bornes
    bornes = []
    indices = []
    tol_bornes = constantes['tol_bornes']
    below_tol = False

    i = 0
    while i < len(angles):
        if not below_tol:
            if answer[i] <= tol_bornes:
                bornes.append([angles[i], -1])
                print("appended")
                indices.append([i, -1])
                below_tol = True
        if below_tol:
            if answer[i] > tol_bornes:
                bornes[-1][1] = angles[i]
                indices[-1][1] = i
                below_tol = False
        i += 1

    if bornes == []:
        return -1
    if bornes[-1][1] == -1:
        bornes[-1][1] = 2*np.pi
        indices[-1][1] = constantes['nbAngles']

    print(bornes)
    print(indices)
    for i in bornes:
        print(str(i) + " : " + str(minimumDistancePourAngleVaisseau(i[0], angleAst, derivee, angles)) + " " + str(minimumDistancePourAngleVaisseau(i[1], angleAst, derivee, angles)))

    # plt.axis('equal')

    # plt.plot(angles, answer, 'g')

    # plt.title("Space simulation")
    # plt.xlabel("Angle")
    # plt.ylabel("minimumDistanceVaisseau")

    # plt.show()

    state = 2
    max_i = len(bornes)
    i = 0
    while state != 0 and i < max_i:
        reponse, state = secante(lambda angleVaisseau: minimumDistancePourAngleVaisseau(angleVaisseau, angleAst, derivee, angles), bornes[i][0], bornes[i][1], constantes['tol'])
        print("State " + str(i) + " : " + str(state))
        i += 1

    constantes = cte

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