import matplotlib.pyplot as plt
import numpy as np

# Funções definidas pelos criadores

def Matriz_A(ml, nl):
    h = 1/ml
    k = 0.5/nl
    Lambda = k/h**2
    MK = np.zeros((ml-1, ml-1))
    for i in range(0, ml-2):
        MK[i][i] = 1 - 2*Lambda
        MK[i][i + 1] = Lambda
        MK[i + 1][i] = Lambda
    MK[ml - 2][ml - 2] = 1 - 2*Lambda
    return MK

def Vetor_nos(ml):
    vnos = np.linspace(1/ml, 1-1/ml, ml-1)
    return vnos

def Vetor_condicaoinicial(vetorl):
    for i in range(0, len(vetorl)):
        vetorl[i] = np.sin(np.pi*vetorl[i])
    return vetorl

# Equação da Difusão do Calor -->> Ut = kUxx <<-- K=1

m = int(input('Digite o número de intervalos (espaço) M = '))
n = int(input('Digite o número de intervalos (tempo) N = '))

# Programa Principal %%% Euler Progressivo %%%%

matriz_A = Matriz_A(m, n)

vetor_nos = Vetor_nos(m)

vetor_ci = Vetor_condicaoinicial(vetor_nos)

for i in range(0, n):
    vetor_ci = np.dot(matriz_A, vetor_ci)
#print(vetor_ci)

# Gráfico
z = np.linspace(0, 1, 100)
plt.plot(z, np.exp(-(np.pi**2)*0.5)*np.sin(np.pi*z))

plt.scatter(Vetor_nos(m), vetor_ci, edgecolors='r')
plt.show()
