from numpy.linalg import det


def linearDependence2(v1, v2):
    flag = True
    r = []
    for i in range(0, 3):
        if v2[i] != 0:
            r.append(v1[i]/v2[i])
    if len(r) > 1:
        for i in range(0, len(r)-1):
            for j in range(i+1, len(r)):
                if abs(r[i]-r[j]) > 1/1000:
                    flag = False
                    break
    
    for i in range(0, 3):
        if (v2[i] == 0) and (v1[i] != 0):
            flag = False
            
    return flag


def linearDependence3(v1, v2, v3):
    M = [v1, v2, v3]
    if det(M) == 0:
        return True
    else:
        return False


def checkE1(V):
    anyTwo = True
    for i in range(0, 2):
        for j in range(i+1, 3):
            if linearDependence2(V[i], V[j]):
                anyTwo = False
                break
    anyThree = linearDependence3(V[0], V[1], V[2])
    return (anyTwo and anyThree)
