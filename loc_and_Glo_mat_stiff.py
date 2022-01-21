import numpy as np
from matplotlib import pyplot as plt


def linear_func(L, n, DBC, xsigm, sigm, fx, E1, E2, A, N):
    l = L/n #Расчет длины элементов
    local_matrix_1 = np.array([
        [A*E1/l,-A*E1/l],
        [-A*E1/l,A*E1/l]
    ]) #Задаем локальную матрицу жесткости для 1 части стержня
    local_matrix_2 = np.array([
        [A*E2/l,-A*E2/l],
        [-A*E2/l,A*E2/l]
    ]) #Задаем локальную матрицу жесткости для 2 части стержня
    Q=np.zeros((n+1,n+1)) #Задаем нулевую глобальную матрицу (все элементы матрицы равны 0)
    #-------------------------Заполнение глобальной матрицы-----------------------------------#
    shift = [0,0]
    for k in range(n//2):
      for i in range(2):
          for j in range(2):
              Q[shift[0]+i][shift[1]+j] += local_matrix_1[i][j]
      for i in range(2):
          shift[i] += 1
    for k in range(n//2, n):
      for i in range(2):
          for j in range(2):
              Q[shift[0]+i][shift[1]+j] += local_matrix_2[i][j]
      for i in range(2):
          shift[i] += 1
    print(Q)
    #-------------------------Матрица сил----------------------------------------#
    F = np.zeros((n+1))
    for i in range(n+1):
        for j in range(len(xsigm)):
            if i==(xsigm[j]):
                F[i-1]+=fx*l+sigm[j]
    for k in range(n+1):
        if F[k]==0:
            F[k]+=fx*l
    F = F.reshape(-1,1)
    #-------------------------Матрица перемещений--------------------------------#
    U=(np.linalg.pinv(Q)).dot(F)
    for i in range(len(U)):
        if i==DBC:
            U[i-1]=0
    print(U)
    #-------------------------График (sigma,x) (displacement,x)--------------------------------#
    for i in range(n+1):
        N.append(i)
    P = np.zeros(n+1)
    for i in range(n+1):
        for j in range(len(xsigm)):
            if i == (xsigm[j]):
                P[i-1]+=sigm[j]
    fig, axes = plt.subplots()
    axes.plot(N,U)
    axes.plot(N,P)
    plt.show()
    print('\n\n\n')


def quadratic_func(L, n, DBC, xsigm, sigm, fx, E1, E2, A, N):
    l = L / n  # Длины элементов
    local_matrix_1 = np.array([
        [A * E1 * 7 / (6 * l), -A * E1 * 4 / (3 * l), A * E1 / (6 * l)],
        [-A * E1 * 4 / (3 * l), A * E1 * 8 / (3 * l), -A * E1 * 4 / (3 * l)],
        [A * E1 / (6 * l), -A * E1 * 4 / (3 * l), A * E1 * 7 / (6 * l)]
    ])
    local_matrix_2 = np.array([
        [A * E2 * 7 / (6 * l), -A * E2 * 4 / (3 * l), A * E2 / (6 * l)],
        [-A * E2 * 4 / (3 * l), A * E2 * 8 / (3 * l), -A * E2 * 4 / (3 * l)],
        [A * E2 / (6 * l), -A * E2 * 4 / (3 * l), A * E2 * 7 / (6 * l)]
    ])
    Q = np.zeros((n + 1, n + 1))  # Глобальная матрица со всеми нулями
    # -------------------------Заполнение матрицы-----------------------------------#
    shift = [0, 0]
    for k in range(n // 2):
        for i in range(2):
            for j in range(2):
                Q[shift[0] + i][shift[1] + j] += local_matrix_1[i][j]
        for i in range(2):
            shift[i] += 1
    for k in range(n // 2, n):
        for i in range(2):
            for j in range(2):
                Q[shift[0] + i][shift[1] + j] += local_matrix_2[i][j]
        for i in range(2):
            shift[i] += 1
    print(Q)
    # -------------------------Матрица сил----------------------------------------#
    F = np.zeros((n + 1))
    for i in range(n + 1):
        for j in range(len(xsigm)):
            if i == (xsigm[j]):
                F[i - 1] += fx * l * l + sigm[j]
    for k in range(n + 1):
        if F[k] == 0:
            F[k] += fx * l * l
    F = F.reshape(-1, 1)
    # -------------------------Матрица перемещений--------------------------------#
    U = (np.linalg.pinv(Q)).dot(F)
    for i in range(len(U)):
        if i == DBC:
            U[i] = 0
    print(U)
    # -------------------------График (sigma,x) (displacement,x)--------------------------------#
    for i in range(n + 1):
        N.append(i)
    print(N)
    P = np.zeros(n + 1)
    for i in range(n + 1):
        for j in range(len(xsigm)):
            if i == (xsigm[j]):
                P[i - 1] += sigm[j]
    fig, axes = plt.subplots()
    axes.plot(N, U)
    axes.plot(N, P)
    plt.show()
    print('\n\n\n')


def start(variant):
    print('Начальные значения:\n')
    xsigm = []
    sigm = []
    N = []
    L = int(input('Длина балки: '))
    n = int(input('Количество элементов: '))
    DBC = int(input('Номер элемента, равный 0: '))
    nsigma = int(input('Количество элементов с сигмой: '))
    #
    for i in range(nsigma):
        a = int(input('Номер элементов с сигмой: '))
        xsigm.append(a)
    for i in range(nsigma):
        v = int(input(f'Сигма на {xsigm[i]} элементе равна: '))
        sigm.append(v)
    fx = int(input('Массовая силу: '))
    E1 = float(input('Модуль Юнга в первой части балки: '))
    E2 = float(input('Модуль Юнга во второй части балки: '))
    A = int(input('Площадь: '))
    if (variant == 1):
        linear_func(L, n, DBC, xsigm, sigm, fx, E1, E2, A, N)
    elif (variant == 2):
        quadratic_func(L, n, DBC, xsigm, sigm, fx, E1, E2, A, N)


if __name__ == "__main__":
    while True:
        variant = int(input('Выберите вариант: 1 - линейная, 2 - квадратная, 3 - выход из программы: '))
        if variant == 1 or variant == 2:
            start(variant)
        else:
            break