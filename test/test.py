from Trajectoire import angle, minimumDistance
from vaisseau import angleVaisseau, minimumDistanceVaisseau
import numpy as np

from valeurs import vitesse, constantes

angleAst = 5.106550164787945

test = angleVaisseau([0, -1.451*100], [vitesse(3.023*1e4), 0], [100, 200], constantes['vitesseAsteroide'], angleAst, constantes['vitesseVaisseau'])
print(test)
dist = minimumDistanceVaisseau(angleAst, test)
print(dist)