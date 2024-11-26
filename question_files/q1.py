class ComplexNumber:
    def __init__(self, real, imag):
        self.real = real
        self.imag = imag

    def __add__(self, other):
        return ComplexNumber(self.real + other.real, self.imag + other.imag)

    def __mul__(self, other):
        return ComplexNumber(
            self.real * other.real - self.imag * other.imag,
            self.real * other.imag + self.imag * other.real,
        )

    def __truediv__(self, other):
        if other.real == 0 and other.imag == 0:
            raise ZeroDivisionError("Division by zero is not allowed.")
        conjugate = other.conjugate()
        numerator = self * conjugate
        denominator = other.real**2 + other.imag**2
        return ComplexNumber(numerator.real / denominator, numerator.imag / denominator)

    def __abs__(self):
        return (self.real**2 + self.imag**2) ** 0.5

    def conjugate(self):
        return ComplexNumber(self.real, -self.imag)

    def __str__(self):
        return f"{self.real} + {self.imag}i"


class Vector:
    
    def __init__(self, field, length):
        self.field = field  
        self.length = length
        self.values = [0] * length

    def set_values(self, values):
        if len(values) != self.length:
            raise ValueError("Input values do not match vector length.")
        self.values = values

    def __add__(self, other):
        if self.length != other.length:
            raise ValueError("Vectors must have the same length for addition.")
        return Vector(self.field, self.length, [self.values[i] + other.values[i] for i in range(self.length)])

    def __str__(self):
        return f"Vector({self.values})"
    
class Matrix:

    def __init__(self, field, rows, cols):
        if field not in ["real", "complex"]:
            raise ValueError("Field must be either 'real' or 'complex'.")
        if rows <= 0 or cols <= 0:
            raise ValueError("Matrix dimensions must be positive integers.")
        self.field = field
        self.rows = rows
        self.cols = cols
        self.entries = [[0] * cols for _ in range(rows)]  

    def set_entries(self, values):
        if len(values) != self.rows or any(len(row) != self.cols for row in values):
            raise ValueError("Input dimensions do not match matrix dimensions.")
        if self.field == "real" and not all(isinstance(v, (int, float)) for row in values for v in row):
            raise ValueError("All entries must be real numbers.")
        if self.field == "complex" and not all(isinstance(v, ComplexNumber) for row in values for v in row):
            raise ValueError("All entries must be complex numbers.")
        self.entries = values

    @classmethod
    def from_columns(cls, field, vectors):
        if not all(v.field == field for v in vectors):
            raise ValueError("All vectors must belong to the specified field.")
        if not all(v.length == vectors[0].length for v in vectors):
            raise ValueError("All vectors must have the same length.")
        rows = vectors[0].length
        cols = len(vectors)
        matrix = cls(field, rows, cols)
        matrix.set_entries([[v.values[i] for v in vectors] for i in range(rows)])
        return matrix

    def __add__(self, other):
        if self.rows != other.rows or self.cols != other.cols:
            raise ValueError("Matrices must have the same dimensions for addition.")
        return Matrix(self.field, self.rows, self.cols, [
            [self.entries[i][j] + other.entries[i][j] for j in range(self.cols)] for i in range(self.rows)
        ])

    def __mul__(self, other):
        if self.cols != other.rows:
            raise ValueError("Matrix multiplication not possible with these dimensions.")
        result = Matrix(self.field, self.rows, other.cols)
        result.set_entries([
            [sum(self.entries[i][k] * other.entries[k][j] for k in range(self.cols)) for j in range(other.cols)]
            for i in range(self.rows)
        ])
        return result

    def get_row(self, index):
        if index < 0 or index >= self.rows:
            raise ValueError("Row index out of bounds.")
        return self.entries[index]

    def get_column(self, index):
        if index < 0 or index >= self.cols:
            raise ValueError("Column index out of bounds.")
        return [self.entries[i][index] for i in range(self.rows)]

    def transpose(self):
        transposed = Matrix(self.field, self.cols, self.rows)
        transposed.set_entries([[self.entries[j][i] for j in range(self.rows)] for i in range(self.cols)])
        return transposed

    def conjugate(self):
        if self.field != "complex":
            raise ValueError("Conjugate is only defined for complex matrices.")
        conjugated = Matrix(self.field, self.rows, self.cols)
        conjugated.set_entries([
            [entry.conjugate() for entry in row] for row in self.entries
        ])
        return conjugated

    def transpose_conjugate(self):
        return self.transpose().conjugate()

    def __str__(self):
        return "\n".join([" ".join(map(str, row)) for row in self.entries])

