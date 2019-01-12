import json


def read_settings(path):
    with open(path) as f:
        data = json.load(f)
    return data
