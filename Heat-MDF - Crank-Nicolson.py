"""Método de Crank-Nicolson"""

import matplotlib.pyplot as plt
import numpy as np

# Funções definidas pelos criadores

def Matriz_A(ml, nl):
    h = 1/ml
    k = 0.5/nl
    Lambda = k/h**2 # Foi considerado alpha = 1
    MK = np.zeros((ml-1, ml-1))
    for i in range(0, ml-2):
        MK[i][i] = 1 + Lambda
        MK[i][i + 1] = - 0.5*Lambda
        MK[i + 1][i] = - 0.5*Lambda
    MK[ml - 2][ml - 2] = 1 + Lambda
    return MK

def Matriz_B(ml, nl):
    h = 1/ml
    k = 0.5/nl
    Lambda = k/h**2 # Foi considerado alpha = 1
    MK = np.zeros((ml-1, ml-1))
    for i in range(0, ml-2):
        MK[i][i] = 1 - Lambda
        MK[i][i + 1] = 0.5*Lambda
        MK[i + 1][i] = 0.5*Lambda
    MK[ml - 2][ml - 2] = 1 - Lambda
    return MK

def Vetor_nos(ml):
    vnos = np.linspace(1/ml, 1-1/ml, ml-1)
    return vnos

def Vetor_condicaoinicial(vetorl):
    for i in range(0, len(vetorl)):
        vetorl[i] = np.sin(np.pi*vetorl[i])
    return vetorl

# Entrada de Dados para a equação do Calor Transiente -->> Ut = kUxx <<-- Supondo K=1

m = int(input('Digite a quantidade de intervalos (espaço) m = '))
n = int(input('Digite a quantidade de intervalos (tempo) n = '))

# Programa Principal %%% Método de Crank-Nicolson %%%%

matriz_A = Matriz_A(m, n)

matriz_B = Matriz_B(m, n)

vetor_nos = Vetor_nos(m)

vetor_ci = Vetor_condicaoinicial(vetor_nos)

for i in range(0, n):
    vetor_ci = np.dot(matriz_B, vetor_ci)
    vetor_ci = np.linalg.solve(matriz_A, vetor_ci)
#print(vetor_ci)

# Gráfico
z = np.linspace(0, 1, 100)
plt.plot(z, np.exp(-(np.pi**2)*0.5)*np.sin(np.pi*z))

plt.scatter(Vetor_nos(m), vetor_ci, edgecolors='r')
plt.show()