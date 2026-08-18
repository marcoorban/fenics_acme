import yaml

meshFile = "mesh.yaml"
paramsFile = "params.yaml"

with (
    open(meshFile, newline="", mode="r") as m,
    open(paramsFile, newline="", mode="r") as p,
):
    meshData = yaml.load(m, Loader=yaml.Loader)
    params = yaml.load(p, Loader=yaml.Loader)
    print(params)
