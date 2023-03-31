# Calcul de la constante de Cavendish de façon précise
def cavendish():
    si = 6.6743*1e-11
    num = ((365*24*3600)**2)*1e24
    den = 1e27
    frac = num/den
    return si*frac

# Transforme une vitesse SI en vitesse astronomique
def vitesse(vit):
    num = 365*24*3600
    den = 1e9
    frac = num/den
    return vit*frac

constantes = {
    'nbMaxIter': 2**8,
    'masseTerre': 5.972,
    'masseAsteroide': 1e-14,
    'masseSoleil': 1.988e6,
    'cavendish': cavendish(),
    'pas': 1, # Pas minimal pour euler ou solve_ivp
    'tol': 1e-1, # Tolérance pour bissection et sécante
    'tol_bornes': 1000, # Tolérance pour trouver les bornes dans la fonction angle
    'r-tol': 1e-2, # rtol et atol de solve_ivp
    'intervalleSimulation': [0, 1],
    'valeursInit': [
        [0, -1.451*100, vitesse(3.023*1e4), 0, 5.972], # Terre
        [100, 200, 0, 0, 1e-14], # Astéroïde
    ],
    'vitesseAsteroide': vitesse(4*1e4),
    'vitesseVaisseau': vitesse(12000),
    'masseVaisseau': 1e-20
}