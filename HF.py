import numpy as np
import math


pi = 4*math.atan(1)


def overlap(basis):
    size = len(basis)
    S = []
    for i in range(size):
        element = []
        for j in range(size):
           result = (pi/(basis[i][0]+basis[j][0]))**1.5
           element.append(result)
        S.append(element)
    return S


def interaction(basis):
    size = len(basis)
    Q = []
    for i in range(size):
        e1 = []
        for j in range(size):
            e2 = []
            for k in range(size):
                e3 = []
                for l in range(size):
                    result = 2*pi**2.5/((basis[i][0]+basis[k][0])*(basis[j][0]+basis[l][0])*math.sqrt(basis[i][0]+basis[j][0]+basis[k][0]+basis[l][0]))
                    e3.append(result)
                e2.append(e3)
            e1.append(e2)
        Q.append(e1)
    return Q


def simpleh(basis,znucl):
    size = len(basis)
    h = []
    for i in range(size):
        element = []
        for j in range(size):
            t = 3*basis[i][0]*basis[j][0]*pi**1.5/(basis[i][0]+basis[j][0])**2.5
            a = -2*znucl*pi/(basis[i][0]+basis[j][0])
            result = t+a
            element.append(result)
        h.append(element)
    return h


def scf(basis, znucl, K):
    size = len(basis)
    h = simpleh(basis, znucl)
    Q = interaction(basis)
    S = overlap(basis)
    S0 = np.linalg.inv(np.array(S))
    C = []
    for k in range(K):
        C.append([1]*size)
    Eg_new = 0
    dEg = 1
    while dEg>0.000001:
        energy1 = 0
        energy2 = 0
        # 归一化C
        for k in range(K):
            result = 0
            for i in range(size):
                for j in range(size):
                    result = result+C[k][i]*C[k][j]*S[i][j]
            C[k] = [x/math.sqrt(result) for x in C[k]]
        # 建立密度矩阵P
        P = []
        for i in range(size):
            element = []
            for j in range(size):
                result = 0
                for k in range(K):
                    result = result+2*C[k][i]*C[k][j]
                element.append(result)
            P.append(element)
        # 建立Fock矩阵并计算总能量
        F = []
        for i in range(size):
            element = []
            for j in range(size):
                result = 0
                energy1 = energy1+P[i][j]*h[i][j]
                for k in range(size):
                    for l in range(size):
                        result = result+P[l][k]*(Q[i][k][j][l]-Q[i][k][l][j]/2)
                        energy2 = energy2+0.5*P[i][j]*P[l][k]*(Q[i][k][j][l]-Q[i][k][l][j]/2)
                element.append(h[i][j]+result)
            F.append(element)
        Eg = Eg_new
        Eg_new = energy1+energy2
        dEg = abs(Eg-Eg_new)
        # 求解本征方程将前k个能量最低的本征值对应的本征矢量作为下一步迭代的C
        F0 = np.array(F)
        C0 = np.linalg.eig(np.matmul(S0, F0))
        C2 = C0[0].tolist()
        C3 = C0[0].tolist()
        C2.sort()
        for k in range(K):
            ind = C3.index(C2[k])
            for i in range(size):
                C[k][i] = C0[1][i][ind]
        print(Eg_new)
    return C0


if __name__ == '__main__':
    basis = [[0.298073, 0, 0], [1.242567, 0, 0], [5.782948, 0, 0], [38.474970, 0, 0]]
    C = scf(basis,znucl = 2,K = 1)
