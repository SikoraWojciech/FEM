from helper.config import read_settings
from controller.ModelsManager import ModelsManager
from shape_functions.matrices import H_matrix_local, C_matrix, H_BC_matrix

def main():
    settings = read_settings()
    models_manager = ModelsManager(settings)
    models_manager.create_grid()
    # models_manager.grid.print()
    x_lst = [0, 0.025, 0.025, 0]
    y_lst = [0, 0, 0.025, 0.025]
    # print(C_matrix(x_lst, y_lst, settings["c"], settings["ro"]))
    print(H_BC_matrix([1, 3], x_lst, y_lst, settings["alfa"]))


if __name__ == '__main__':
    main()
