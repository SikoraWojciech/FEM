from helper.config import read_settings
from controller.ModelsManager import ModelsManager
import sys


def main():

    if len(sys.argv) < 3:
        raise Exception('Not enough params')

    settings = read_settings(sys.argv[1])
    settings["img_path"] = sys.argv[2]

    models_manager = ModelsManager(settings)
    models_manager.create_grid()
    models_manager.grid.heat(settings)


if __name__ == '__main__':
    main()
