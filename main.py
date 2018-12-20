from helper.config import read_settings
from controller.ModelsManager import ModelsManager


def heat():
    return


def main():
    settings = read_settings()
    models_manager = ModelsManager(settings)
    models_manager.create_grid()
    heat()


if __name__ == '__main__':
    main()
