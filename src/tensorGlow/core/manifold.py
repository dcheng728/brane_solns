import sympy as sp


class Operator():
    def __init__(self, name, manifold):
        self.name = name
        self.manifold = manifold

    def __str__(self):
        return f"{self.name}"

class Manifold:
    def __init__(self, signature):
        self.num_negative_dim = signature[0]
        self.num_positive_dim = signature[1]
        self.coordinates = [Operator(f"x{i}", self) for i in range(self.num_negative_dim + self.num_positive_dim)]

    def get_coordinates(self):
        return self.coordinates

    def __str__(self):
        return f"Manifold with {self.num_negative_dim} negative dimensions and {self.num_positive_dim} positive dimensions"
    

if __name__ == "__main__":
    m = Manifold((1, 3))
    print(m)

    print("Coordinates:", m.get_coordinates())