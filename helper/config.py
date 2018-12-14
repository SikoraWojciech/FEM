import json


def read_settings():
    with open('./settings.json') as f:
        data = json.load(f)
    return data
