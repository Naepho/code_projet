import numpy as np

# Passe de coordonnées polaires aux cartésiennes
def compVect(angle, norme):
    compX = np.cos(angle)*norme
    compY= np.sin(angle)*norme
    
    return [compX,compY]