import numpy as np
from valeurs import constantes

def bissection(f, x0, x1, tol):
    # print("Début de la bissection")

    try:
        fx0 = f(x0)
        fx1 = f(x1)
    except:
        print("Une erreur est survenue, la fonction n'est porbablement pas continue")
        return [0, -1]

    stock = 0

    # Vérification des hypothèses sur les valeurs d'entrée
    if not ((fx0 > 0 and fx1 < 0) or (fx0 < 0 and fx1 > 0)):
        print("Les absices rentrées ne satisfont pas aux hypothèses de la bissection")
        x = 0
        statut = 1
        return [x, statut]

    # Vérification sur f(x0) et f(x1) pour voir si l'un des deux n'est pas la réponse
    if np.abs(fx0) < tol or np.abs(fx1) < tol:
        x = x0 if np.abs(fx0) < tol else x1
        statut = 0
        return [x, statut]

    # permutation des absices pour pouvoir appliquer toujours le même procédé
    # x0 est négatif et x1 est positif
    if fx1 < 0:
        stock = x1
        x1 = x0
        x0 = stock

    # Algorithme en lui-même
    absmilieu = (x1 + x0) / 2
    try:
        fmilieu = f(absmilieu)
    except:
        print("Une erreur est survenue, la fonction n'est porbablement pas continue")
        return [0, -1]
    count = 0

    while np.abs(fmilieu) > tol and count < constantes["nbMaxIter"]:
        if fmilieu < 0:
            x0 = absmilieu
        else:
            x1 = absmilieu
        absmilieu = (x1 + x0) / 2
        try:
            fmilieu = f(absmilieu)
        except:
            print(
                "Une erreur est survenue, la fonction n'est porbablement pas continue"
            )
            return [0, -1]
        count += 1

    if count == constantes["nbMaxIter"]:
        return [0, -1]
    return [absmilieu, 0]


def secante(f, x0, x1, tol):
    # print("Début de la sécante")

    try:
        fx0 = f(x0)
        fx1 = f(x1)
    except:
        print("Une erreur est survenue, la fonction n'est porbablement pas continue")
        return [0, -1]
    
    # Vérification de l'hypothèse sur les ordonnées des points 
    if fx0 == fx1:
        print("Ordonnées égales, transmettez de nouvelles valeurs pour les abssices")
        return [0, 1]

    # Vérification sur f(x0) et f(x1) pour voir si l'un des deux n'est pas la réponse
    if np.abs(fx0) < tol or np.abs(fx1) < tol:
        x = x0 if np.abs(fx0) < tol else x1
        statut = 0
        return [x, statut]

    # Algorithme en lui-même

    # Initialisation
    absnew = x1 - (fx1 * (x1 - x0))/(fx1 - fx0)
    try:
        fnew = f(absnew)
    except:
        print("Une erreur est survenue, la fonction n'est porbablement pas continue")
        return [0, -1]
    count = 0

    # x0 est la borne la plus ancienne et x1 est la plus récente
    # Boucle principale de l'algorithme
    while np.abs(fnew) > tol and count < constantes["nbMaxIter"]:
        absnew = x1 - (fx1 * (x1 - x0))/(fx1 - fx0)
        try:
            fnew = f(absnew)
        except:
            print("Une erreur est survenue lors du calcul de la valeur de la fonction.")
            print("Elle n'est surement pas continue")
            return [0, -1]

        if fnew == fx1:
            print("Une nouvelle abscisse est égale à l'autre borne, on ne va pas diviser par zero")
            return [0, -1]
        
        # Permutation pour que x0 reste la borne la plus ancienne
        x0 = x1
        fx0 = fx1
        x1 = absnew
        fx1 = fnew

        count += 1

    if count == constantes["nbMaxIter"]:
        return [0, -1]
    
    return [absnew, 0]

# Si la racine existe, on la trouvera
def recursive(f, x0, x1, tol, tol_bornes, evolTol, nbVal):
    bornes = []
    indices = []
    below_tol = False
    pas = (x1-x0)/nbVal

    answer = np.zeros(nbVal)
    for i in range(nbVal):
        answer[i] = f(x0 + i*pas)
        if np.abs(answer[i]) <= tol:
            return x0+i*pas, 0

    i = 0
    while i < nbVal:
        if not below_tol:
            if np.abs(answer[i]) <= tol_bornes:
                bornes.append([x0 + (i-1)*pas, -1])
                if i == 0:
                    bornes[-1][0] = x0
                indices.append([i, -1])
                below_tol = True
        if below_tol:
            if np.abs(answer[i]) > tol_bornes:
                bornes[-1][1] = x0 + i*pas
                indices[-1][1] = i
                below_tol = False
        i += 1

    if bornes == []:
        return 0, -1
    if bornes[-1][1] == -1:
        bornes[-1][1] = x1
        indices[-1][1] = nbVal
    
    for i in range(len(bornes)):
        reponse_potentielle, state = recursive(f, bornes[i][0], bornes[i][1], tol, tol_bornes*evolTol, evolTol, nbVal)
        if state != -1:
            return reponse_potentielle, 0
    return 0, -1
