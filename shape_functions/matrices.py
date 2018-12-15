from models.models import PointKsiEta
from shape_functions.functions import *
from math import *
from numpy import *


def dN_dksi_matix(p):
    result = []
    for i in range(4):
        result.append(dN_dksi(i, p.eta))
    return result


def dN_deta_matix(p):
    result = []
    for i in range(4):
        result.append(dN_deta(i, p.eta))
    return result


def Jacobian_matrix(x_lst, y_lst, dN_dksi_lst, dN_deta_lst):
    dx_dksi = 0
    dx_deta = 0
    dy_dksi = 0
    dy_deta = 0
    # dx/dksi = dN1/dksi * x1 + .. + dN4/dksi * x4
    # dy/dksi = dN1/deta * x1 + .. + dN4/deta * x4
    # analogicznie dla y
    for i in range(4):
        dx_dksi += dN_dksi_lst[i] * x_lst[i]
        dx_deta += dN_deta_lst[i] * x_lst[i]
        dy_dksi += dN_dksi_lst[i] * y_lst[i]
        dy_deta += dN_deta_lst[i] * y_lst[i]
    matrix_tmp = matrix([[dx_dksi, dy_dksi],
                         [dx_deta, dy_deta]])
    return matrix_tmp


def dN_dx_matrix(inverse_jacobian, dN_dksi_lst, dN_deta_lst):
    # Operujemy na gornych pozycjach odwroconej macierzy Jacobianu
    # Wzor na dN1/dx = J-1(x, ksi) * dN1/dksi + J-1(x, eta) * dN1/deta
    result = zeros([4])
    for i in range(4):
        result.itemset(i, inverse_jacobian.item(0) * dN_dksi_lst[i] + inverse_jacobian.item(1) * dN_deta_lst[i])
    return result


def dN_dy_matrix(inverse_jacobian, dN_dksi_lst, dN_deta_lst):
    # Operujemy na dolnych pozycjach odwroconej macierzy Jacobianu
    # Wzor na dN1/dx = J-1(x, ksi) * dN1/dksi + J-1(x, eta) * dN1/deta
    result = zeros([4])
    for i in range(4):
        result.itemset(i, inverse_jacobian.item(2) * dN_dksi_lst[i] + inverse_jacobian.item(3) * dN_deta_lst[i])
    return result


def H_matrix_local(x_lst, y_lst, k):
    # 4 punkty calkowania - uklad ksi i eta
    p = [
        PointKsiEta(-1 / sqrt(3), -1 / sqrt(3)),
        PointKsiEta(1 / sqrt(3), -1 / sqrt(3)),
        PointKsiEta(1 / sqrt(3), 1 / sqrt(3)),
        PointKsiEta(-1 / sqrt(3), 1 / sqrt(3))
    ]

    matrix_multiplied_by_k = []
    for i in range(4):
        # Wyznaczamy listy 4x1 dN/dksi i dN/deta dla kazdego punktu calkowania
        dN_dksi_lst = dN_dksi_matix(p[i])
        dN_deta_lst = dN_deta_matix(p[i])

        # Wykorzystujemy je do wyliczenia Jacobianu
        jacobian = Jacobian_matrix(x_lst, y_lst, dN_dksi_lst, dN_deta_lst)

        # Odwracamy Jacobian by wyznaczyc dN/dx i dN/dy
        inverse_jacobian = linalg.inv(jacobian)

        # Wyznacznik bedzie nam potrzebny by powiedziec na ile oszacowalismy wynik
        det_jacobian = linalg.det(jacobian)

        # Mnozymy listy dN/dx i (dN/dx)T przez wyznacznik
        dN_dx = dN_dx_matrix(inverse_jacobian, dN_dksi_lst, dN_deta_lst)
        dN_dx_T = dN_dx.reshape(4, 1)
        dN_dx_matrix_multiplied_by_det = (dN_dx * dN_dx_T) * det_jacobian

        # Mnozymy listy dN/dy i (dN/dy)T przez wyznacznik
        dN_dy = dN_dy_matrix(inverse_jacobian, dN_dksi_lst, dN_deta_lst)
        dN_dy_T = dN_dy.reshape(4, 1)
        dN_dy_matrix_multiplied_by_det = (dN_dy * dN_dy_T) * det_jacobian

        # Mnozymy nasze macierze 4x4 razy wspolczynnik przewodzenia k
        matrix_multiplied_by_k.append((dN_dy_matrix_multiplied_by_det + dN_dx_matrix_multiplied_by_det) * k)

    # Wyznaczamy macierz H 4x4 poprzez sumowanie poszczegolnych wartosci z macierzy pomnozonej przez k
    H_matrix = zeros([4, 4])
    for i in range(4):
        for j in range(4):
            val_tmp = 0
            for integration_point in range(4):
                val_tmp += matrix_multiplied_by_k[integration_point].item((i, j))
            H_matrix.itemset((i, j), val_tmp)

    return H_matrix


def C_matrix_local(x_lst, y_lst, c, ro):
    # Macierz C : Calka(c * ro * {N} * {N}T)
    # Mnozymy przez wyznacznik Jakobianu zeby miec wynik calki (powiedziec na ile oklamalismy)
    p = [
        PointKsiEta(-1 / sqrt(3), -1 / sqrt(3)),
        PointKsiEta(1 / sqrt(3), -1 / sqrt(3)),
        PointKsiEta(1 / sqrt(3), 1 / sqrt(3)),
        PointKsiEta(-1 / sqrt(3), 1 / sqrt(3))
    ]

    N_matrices_multiplied = []
    for i in range(4):
        # Wyznaczamy listy 4x1 dN/dksi i dN/deta dla kazdego punktu calkowania
        dN_dksi_lst = dN_dksi_matix(p[i])
        dN_deta_lst = dN_deta_matix(p[i])

        # Wykorzystujemy je do wyliczenia Jacobianu
        jacobian = Jacobian_matrix(x_lst, y_lst, dN_dksi_lst, dN_deta_lst)

        # Wyznacznik bedzie nam potrzebny by powiedziec na ile oszacowalismy wynik
        det_jacobian = linalg.det(jacobian)

        # Obliczanie {N}*{N}T * c * ro * det Jacobian
        N_matrix = zeros(4)
        for j in range(4):
            N_matrix.itemset(j, N(j, p[i].ksi, p[i].eta))
        N_matrix_T = N_matrix.reshape(4, 1)
        N_matrix_multiplied = N_matrix * N_matrix_T

        N_matrices_multiplied.append(N_matrix_multiplied * c * ro * det_jacobian)

    # Obliczanie macierzy C na podstawie punktow calkowania
    # Macierz C jest suma odpowiednich wartosci z macierzy poszczegolnych punktow calkowania
    C_matrix = zeros([4, 4])
    for i in range(4):
        for j in range(4):
            val_tmp = 0
            for integration_point in range(4):
                val_tmp += N_matrices_multiplied[integration_point].item(i, j)
            C_matrix.itemset((i, j), val_tmp)
    return C_matrix
