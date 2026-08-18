import io


class FEM_Solver:
    def __init__(self, paramsFile):
        self.setParams(paramsFile)

    def setParams(self, paramsFile):
        params = io.readParams(paramsFile)
        self.physics = params["physics"]
        self.boundary_conditions = params["boundary_conditions"]
        self.nitsche = params["Nitsche"]
        return

    def print_params(self):
        print(self.physics, self.boundary_conditions, self.nitsche)


strong = FEM_Solver("params.yaml")
strong.print_params()
