import yaml


def readParams(paramsFile):
    with open(paramsFile, newline="", mode="r") as p:
        return yaml.load(p, Loader=yaml.Loader)
