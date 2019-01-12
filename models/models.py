import matplotlib.patches as mpatches
from numpy import *
import matplotlib.pyplot as plt
import seaborn as sb
import imageio


class Node:

    def __init__(self, id, x, y, temp):
        self.id = id
        self.x = x
        self.y = y
        self.temp = temp

    def print(self):
        print("ID: {} ({}; {}) temp: {}".format(self.id, self.x, self.y, self.temp))


class Element:

    def __init__(self, id, node1, node2, node3, node4, surface_indexes, k):
        self.id = id
        self.nodes = [node1, node2, node3, node4]
        self.heated_surfaces_indexes = surface_indexes
        self.k = k
        self.H_matrix = zeros([4, 4])
        self.C_matrix = zeros([4, 4])
        self.P_vector = zeros(4)

    def print(self):
        print("ID: {}; Nodes[{}, {}, {}, {}]; k = {}; Surfaces_heated = {}"
              .format(self.id, self.nodes[0].id, self.nodes[1].id,
                      self.nodes[2].id, self.nodes[3].id, self.k, self.heated_surfaces_indexes))


class Grid:

    def __init__(self, settings):
        self.nodes = []
        self.elements = []
        nodes_count = settings["nH"] * settings["nL"]
        self.H_matrix = zeros((nodes_count, nodes_count))
        self.C_matrix = zeros((nodes_count, nodes_count))
        self.P_vector = zeros(nodes_count)

    def __del__(self):
        for node in self.nodes:
            del node
        for element in self.elements:
            del element

    def print(self):
        for node in self.nodes:
            node.print()
        print()
        for element in self.elements:
            element.print()

    def agregate_matrices(self):
        for element in self.elements:
            for i_local in range(4):
                for j_local in range(4):
                    i_global = element.nodes[i_local].id - 1
                    j_global = element.nodes[j_local].id - 1

                    val_H_local = element.H_matrix.item((i_local, j_local))
                    val_H_global = self.H_matrix.item((i_global, j_global)) + val_H_local
                    self.H_matrix.itemset((i_global, j_global), val_H_global)

                    val_C_local = element.C_matrix.item((i_local, j_local))
                    val_C_global = self.C_matrix.item((i_global, j_global)) + val_C_local
                    self.C_matrix.itemset((i_global, j_global), val_C_global)

    def agregate_vector(self):
        for element in self.elements:
            for i_local in range(4):
                i_global = element.nodes[i_local].id - 1

                val_P_local = element.P_vector.item(i_local)
                val_P_global = self.P_vector.item(i_global) + val_P_local
                self.P_vector.itemset(i_global, val_P_global)

    def heat(self, settings):
        nodes_count = len(self.nodes)
        dT = settings["time_step"]

        left_part_of_equation = self.H_matrix + (self.C_matrix / dT)

        max_temp_list = []
        min_temp_list = []

        for interval in range(0, settings["simulation_time"], dT):

            # Wektor t0 - temperatury kazdego z wezlow
            t0_vector = zeros(nodes_count)
            for i in range(nodes_count):
                t0_vector.itemset(i, self.nodes[i].temp)

            right_part_of_equation = -self.P_vector + (self.C_matrix / dT) @ t0_vector

            t1_vector = linalg.solve(left_part_of_equation, right_part_of_equation)

            for i in range(len(t1_vector)):
                self.nodes[i].temp = t1_vector[i]

            # print('Interval {}'.format(interval + dT))
            # print('=========================================')
            # print('Max: {}'.format(max(t1_vector)))
            # print('Min: {}'.format(min(t1_vector)))
            max_temp_list.append(max(t1_vector))
            min_temp_list.append(min(t1_vector))
            if (settings["simulation_time"]/settings["time_step"]) <= 20:
                self.create_heatmap(interval, settings)
            print()
        # self.create_heatmap_animation(settings)
        self.create_chart(max_temp_list, min_temp_list, 'temp_chart', settings)
        self.save_temps_to_csv(max_temp_list, min_temp_list, settings)

    def create_heatmap(self, interval, settings):
        nodes_heatmap = zeros([settings["nH"], settings["nL"]])
        counter = 0
        for j in range(settings["nL"]):
            for i in range(settings["nH"]):
                nodes_heatmap.itemset((i, j), self.nodes[counter].temp)
                counter += 1
        plt.figure(figsize=(settings["nH"], settings["nL"]))
        plt.gca().invert_yaxis()
        sb.heatmap(nodes_heatmap, 0, settings["ambient_temp"], square=True, cmap="coolwarm", yticklabels='',
                   xticklabels='')
        # sb.heatmap(nodes_heatmap, 0, 1000, cmap="coolwarm", yticklabels='',
        #            xticklabels='', square=True)
        plt.title('Time : {}s'.format(interval + settings["time_step"]))
        plt.savefig('heatmaps/heatmap_{}s.jpg'.format(interval + settings["time_step"]))

    def create_chart(self, max_temp_lst, min_temp_lst, chart_name, settings):
        plt.close('all')
        plt.plot(min_temp_lst)
        plt.plot(max_temp_lst)
        plt.xlabel('Time [s]')
        plt.ylabel('Temperature [C]')
        max_temp_lbl = mpatches.Patch(color='orange', label='Max')
        min_temp_lbl = mpatches.Patch(color='blue', label='Min')
        plt.legend(handles=[max_temp_lbl, min_temp_lbl])
        plt.savefig(settings["img_path"] + '{}.jpg'.format(chart_name))

    def create_heatmap_animation(self, settings):
        images = []
        filenames = []
        for interval in range(0, settings["simulation_time"], settings["time_step"]):
            filenames.append('heatmaps/heatmap_{}s.png'.format(interval + settings["time_step"]))
        for filename in filenames:
            images.append(imageio.imread(filename))
        imageio.mimsave(settings["img_path"] + 'heating.gif', images)

    def save_temps_to_csv(self, max_temp_lst, min_temp_lst, settings):
        arr = array([max_temp_lst, min_temp_lst])
        savetxt(settings["img_path"] + '{}.csv'.format('temperatures'), arr, delimiter=";")


class PointKsiEta:
    def __init__(self, ksi, eta):
        self.ksi = ksi
        self.eta = eta
