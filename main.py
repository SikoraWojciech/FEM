from helper.config import read_settings
from controller.ModelsManager import ModelsManager


def main():
    settings = read_settings()
    models_manager = ModelsManager(settings)
    models_manager.create_grid()
    models_manager.grid.heat(settings)


if __name__ == '__main__':
    main()
