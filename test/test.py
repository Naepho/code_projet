from Trajectoire import angle, minimumDistance
from vaisseau import angleVaisseau, minimumDistanceVaisseau
from Calcul_trajectoire import simulation
import numpy as np
import time

from valeurs import vitesse, constantes

# angleAst = 5.106550164787945
angleAst =  74/360*2*np.pi + np.pi

# angleAst = angle([-1.451*100, 0], [0, vitesse(3.023*1e4)], [100, 200], constantes['vitesseAsteroide'])
# print(angleAst)

start = time.time()
test = angleVaisseau([1.451*100, 0], [0, - vitesse(3.023*1e4)], [- 100,  100], constantes['vitesseAsteroide'], angleAst, constantes['vitesseVaisseau'])
print("Time of angle : " + str(time.time() - start))
print(test)
start = time.time()
dist = minimumDistanceVaisseau(angleAst, test)
print(dist)
print("Time of minimum distance : " + str(time.time() - start))
simulation(angleAst, test)
