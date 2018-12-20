from models.models import PointKsiEta, Element
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


def Jacobian_matrix(element, dN_dksi_lst, dN_deta_lst):
    dx_dksi = 0
    dx_deta = 0
    dy_dksi = 0
    dy_deta = 0
    # dx/dksi = dN1/dksi * x1 + .. + dN4/dksi * x4
    # dy/dksi = dN1/deta * x1 + .. + dN4/deta * x4
    # analogicznie dla y
    for i in range(4):
        dx_dksi += dN_dksi_lst[i] * element.nodes[i].x
        dx_deta += dN_deta_lst[i] * element.nodes[i].x
        dy_dksi += dN_dksi_lst[i] * element.nodes[i].y
        dy_deta += dN_deta_lst[i] * element.nodes[i].y
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


def H_matrix(element):
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
        jacobian = Jacobian_matrix(element, dN_dksi_lst, dN_deta_lst)

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
        matrix_multiplied_by_k.append((dN_dy_matrix_multiplied_by_det + dN_dx_matrix_multiplied_by_det) * 30)

    # Wyznaczamy macierz H 4x4 poprzez sumowanie poszczegolnych wartosci z macierzy pomnozonej przez k
    result = zeros([4, 4])
    for i in range(4):
        for j in range(4):
            val_tmp = 0
            for integration_point in range(4):
                val_tmp += matrix_multiplied_by_k[integration_point].item((i, j))
                result.itemset((i, j), val_tmp)
    return result


def C_matrix(element, c, ro):
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
        jacobian = Jacobian_matrix(element, dN_dksi_lst, dN_deta_lst)

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
    # Macierz C jest suma odpowiednich wartosci z macierzy pomnozonej przez transponowana dla punktow calkowania
    result = zeros([4, 4])
    for i in range(4):
        for j in range(4):
            val_tmp = 0
            for integration_point in range(4):
                val_tmp += N_matrices_multiplied[integration_point].item(i, j)
                result.itemset((i, j), val_tmp)
    return result


def H_BC_matrix(element, alfa):
    # Macierz lokalna H dla warunkow brzegowych
    # o - punkty calkowania dla powierzchni
    #
    #           SI = 3
    #         __o_____o__
    #        |           |
    #        o           o
    # SI = 4 |           | SI = 2
    #        o           o
    #        |__o_____o__|
    #
    #           SI = 1

    p = [PointKsiEta(-1 / sqrt(3), -1),
         PointKsiEta(1 / sqrt(3), -1),
         PointKsiEta(1, -1 / sqrt(3)),
         PointKsiEta(1, 1 / sqrt(3)),
         PointKsiEta(1 / sqrt(3), 1),
         PointKsiEta(-1 / sqrt(3), 1),
         PointKsiEta(-1, 1 / sqrt(3)),
         PointKsiEta(-1, -1 / sqrt(3))]

    # Dla kazdego boku mamy 1D wiec wyznacznik macierzy Jakobiego = dlugosc_boku/2
    detJx = element.nodes[1].x - element.nodes[0].x
    detJy = element.nodes[2].y - element.nodes[1].y

    surfaces = [{"SI": 1, "detJ": detJx, "p_range": range(0, 2), "is_heated": 0},
                {"SI": 2, "detJ": detJy, "p_range": range(2, 4), "is_heated": 0},
                {"SI": 3, "detJ": detJx, "p_range": range(4, 6), "is_heated": 0},
                {"SI": 4, "detJ": detJy, "p_range": range(6, 8), "is_heated": 0}]

    for surface_index in element.heated_surfaces_indexes:
        if surface_index == 1:  # SI = 1
            surfaces[0]["is_heated"] = 1
        if surface_index == 2:  # SI = 2
            surfaces[1]["is_heated"] = 1
        if surface_index == 3:  # SI = 3
            surfaces[2]["is_heated"] = 1
        if surface_index == 4:  # SI = 4
            surfaces[3]["is_heated"] = 1

    # Obliczanie {N}*{N}T * alfa dla kazdej powierzchni -> kazdego punktu calkowania
    N_matrices_multiplied = []
    for surface in surfaces:
        N_matrix_for_sum = zeros((4, 4))
        for integral_point_no in surface["p_range"]:
            N_matrix = zeros(4)
            for shape_func_no in range(4):
                N_matrix.itemset(shape_func_no, N(shape_func_no, p[integral_point_no].ksi, p[integral_point_no].eta))
            N_matrix_T = N_matrix.reshape(4, 1)
            N_matrix_for_sum += N_matrix * N_matrix_T * alfa
        N_matrix_for_sum = N_matrix_for_sum * 0.5 * surface["detJ"] * surface["is_heated"]
        N_matrices_multiplied.append(N_matrix_for_sum)

    # Macierz H dla warunkow brzechowych to suma macierzy powstalych z {N}*{N}T * alfa
    result = zeros([4, 4])
    for matrix_index in range(len(N_matrices_multiplied)):
        result += N_matrices_multiplied[matrix_index]
    return result


def H_matrix_local(element, alfa):
    # Macierz lokalna powstaje w wyniku zsumowania macierzy H dla punktow i macierzy H z warunkami brzegowymi
    return H_matrix(element) + H_BC_matrix(element, alfa)


def P_vector(element, ambient_temp, alfa):
    # Sytuacja analogiczna jak w macierzy H z warunkami brzegowymi
    p = [PointKsiEta(-1 / sqrt(3), -1),
         PointKsiEta(1 / sqrt(3), -1),
         PointKsiEta(1, -1 / sqrt(3)),
         PointKsiEta(1, 1 / sqrt(3)),
         PointKsiEta(1 / sqrt(3), 1),
         PointKsiEta(-1 / sqrt(3), 1),
         PointKsiEta(-1, 1 / sqrt(3)),
         PointKsiEta(-1, -1 / sqrt(3))]

    # Dla kazdego boku mamy 1D wiec wyznacznik macierzy Jakobiego = dlugosc_boku/2
    detJx = element.nodes[1].x - element.nodes[0].x
    detJy = element.nodes[2].y - element.nodes[1].y

    surfaces = [{"SI": 1, "detJ": detJx, "p_range": range(0, 2), "is_heated": 0},
                {"SI": 2, "detJ": detJy, "p_range": range(2, 4), "is_heated": 0},
                {"SI": 3, "detJ": detJx, "p_range": range(4, 6), "is_heated": 0},
                {"SI": 4, "detJ": detJy, "p_range": range(6, 8), "is_heated": 0}]

    for surface_index in element.heated_surfaces_indexes:
        if surface_index == 1:  # SI = 1
            surfaces[0]["is_heated"] = 1
        if surface_index == 2:  # SI = 2
            surfaces[1]["is_heated"] = 1
        if surface_index == 3:  # SI = 3
            surfaces[2]["is_heated"] = 1
        if surface_index == 4:  # SI = 4
            surfaces[3]["is_heated"] = 1

    # Obliczanie {N} * -alfa * temp_otoczenia dla kazdej powierzchni -> kazdego punktu calkowania
    P_vectors_multiplied = []
    for surface in surfaces:
        P_vector_for_sum = zeros(4)
        for integral_point_no in surface["p_range"]:
            P_vec = zeros(4)
            for shape_func_no in range(4):
                P_vec.itemset(shape_func_no, N(shape_func_no, p[integral_point_no].ksi, p[integral_point_no].eta))
            P_vector_for_sum += P_vec
        P_vector_for_sum *= 0.5 * surface["detJ"] * surface["is_heated"] * (- alfa) * ambient_temp
        P_vectors_multiplied.append(P_vector_for_sum)

    result = zeros(4)
    for P_vec_index in range(len(P_vectors_multiplied)):
        result += P_vectors_multiplied[P_vec_index]
    return result
