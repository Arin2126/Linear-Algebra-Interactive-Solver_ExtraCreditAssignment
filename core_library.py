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
    
    def append_column(self, column):
        if len(column) != self.rows:
            raise ValueError("The length of the column must match the number of rows in the matrix.")

        augmented_entries = [self.entries[i] + [column[i]] for i in range(self.rows)]
        augmented_matrix = Matrix(self.field, self.rows, self.cols + 1)
        augmented_matrix.set_entries(augmented_entries)
        return augmented_matrix

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





### Core Functions

def compute_rref(matrix):
    """Compute the Reduced Row Echelon Form (RREF) of a matrix."""
    rref_matrix = Matrix("real", matrix.rows, matrix.cols)
    rref_matrix.set_entries([row[:] for row in matrix.entries])  
    pivot_row = 0

    for pivot_col in range(matrix.cols):
        if pivot_row >= matrix.rows:
            break

        max_row = max(range(pivot_row, matrix.rows), key=lambda r: abs(rref_matrix.entries[r][pivot_col]))
        if abs(rref_matrix.entries[max_row][pivot_col]) < 1e-10:
            continue

        rref_matrix.entries[pivot_row], rref_matrix.entries[max_row] = (
            rref_matrix.entries[max_row],
            rref_matrix.entries[pivot_row],
        )

        pivot_value = rref_matrix.entries[pivot_row][pivot_col]
        rref_matrix.entries[pivot_row] = [val / pivot_value for val in rref_matrix.entries[pivot_row]]

        for row in range(matrix.rows):
            if row != pivot_row:
                factor = rref_matrix.entries[row][pivot_col]
                rref_matrix.entries[row] = [
                    rref_matrix.entries[row][i] - factor * rref_matrix.entries[pivot_row][i]
                    for i in range(matrix.cols)
                ]

        pivot_row += 1

    return rref_matrix


def invert_matrix(matrix):
    """Compute the inverse of a square matrix."""
    size = matrix.rows
    augmented = Matrix("real", size, size * 2)
    augmented.set_entries(
        [matrix.entries[i] + [1 if i == j else 0 for j in range(size)] for i in range(size)]
    )
    rref = compute_rref(augmented)
    inverse_entries = [row[size:] for row in rref.entries]
    inverse_matrix = Matrix("real", size, size)
    inverse_matrix.set_entries(inverse_entries)
    return inverse_matrix


def multiply_matrices(matrix1, matrix2):
    """Multiply two matrices."""
    rows = matrix1.rows
    cols = matrix2.cols
    result = Matrix("real", rows, cols)
    result.set_entries(
        [
            [sum(matrix1.entries[i][k] * matrix2.entries[k][j] for k in range(matrix1.cols)) for j in range(cols)]
            for i in range(rows)
        ]
    )
    return result


def multiply_matrix_vector(matrix, vector):
    """Multiply a matrix by a vector."""
    result = Vector("real", matrix.rows)
    result.set_values(
        [sum(matrix.entries[i][j] * vector.values[j] for j in range(matrix.cols)) for i in range(matrix.rows)]
    )
    return result


def print_matrix(matrix):
    """Print a matrix in a readable format."""
    for row in matrix.entries:
        print(" ".join(f"{val:.2f}" for val in row))


def sqrt(x):
    """Approximate square root using Newton's method."""
    if x < 0:
        raise ValueError("Cannot compute square root of a negative number.")
    guess = x / 2
    for _ in range(10):
        guess = (guess + x / guess) / 2
    return guess


def prod(iterable):
    """Compute the product of elements in an iterable."""
    result = 1
    for x in iterable:
        result *= x
    return result


def compute_rref_scaling_factor(matrix):
    """Compute the scaling factor for RREF determinant calculation."""
    size = matrix.rows
    scaling_factor = 1
    for i in range(size):
        for j in range(size):
            if i != j and matrix.entries[i][j] != 0:
                scaling_factor *= matrix.entries[i][j]
    return scaling_factor


def product_of_diagonal_elements(matrix):
    """Compute the product of diagonal elements of a matrix."""
    return prod(matrix.entries[i][i] for i in range(matrix.rows))


def count_permutations(P):
    """Count the number of transpositions in a permutation matrix."""
    visited = [False] * len(P)
    count = 0

    for i in range(len(P)):
        if not visited[i]:
            cycle_length = 0
            j = i
            while not visited[j]:
                visited[j] = True
                j = P[j].index(1)
                cycle_length += 1
            count += cycle_length - 1

    return count

def compute_plu_decomposition(matrix):
    """Compute the PLU decomposition of a square matrix."""
    size = matrix.rows
    if matrix.rows != matrix.cols:
        raise ValueError("PLU decomposition requires a square matrix.")

    P = [[1 if i == j else 0 for j in range(size)] for i in range(size)]
    L = [[0 for _ in range(size)] for _ in range(size)]
    U = [[matrix.entries[i][j] for j in range(size)] for i in range(size)]

    for i in range(size):
        
        max_row = max(range(i, size), key=lambda x: abs(U[x][i]))
        if abs(U[max_row][i]) < 1e-10:
            raise ValueError("Matrix is singular and cannot be decomposed.")

        
        P[i], P[max_row] = P[max_row], P[i]
        U[i], U[max_row] = U[max_row], U[i]

        
        for j in range(i + 1, size):
            factor = U[j][i] / U[i][i]
            L[j][i] = factor
            for k in range(i, size):
                U[j][k] -= factor * U[i][k]

    
    for i in range(size):
        L[i][i] = 1

    
    L_matrix = Matrix("real", size, size)
    U_matrix = Matrix("real", size, size)
    P_matrix = Matrix("real", size, size)
    L_matrix.set_entries(L)
    U_matrix.set_entries(U)
    P_matrix.set_entries(P)

    return P_matrix, L_matrix, U_matrix

def cofactor_expansion(matrix):
    """Compute determinant using the cofactor expansion method."""
    size = matrix.rows
    if matrix.rows != matrix.cols:
        raise ValueError("Cofactor expansion requires a square matrix.")

    
    if size == 1:
        return matrix.entries[0][0]

    
    if size == 2:
        return matrix.entries[0][0] * matrix.entries[1][1] - matrix.entries[0][1] * matrix.entries[1][0]

    
    determinant = 0
    for col in range(size):
        
        sub_matrix = Matrix("real", size - 1, size - 1)
        sub_matrix.set_entries([
            [matrix.entries[i][j] for j in range(size) if j != col]
            for i in range(1, size)
        ])
        cofactor = ((-1) ** col) * matrix.entries[0][col] * cofactor_expansion(sub_matrix)
        determinant += cofactor

    return determinant

def compute_lu_decomposition(matrix):
    """Compute the LU decomposition of a square matrix."""
    size = matrix.rows
    if matrix.rows != matrix.cols:
        raise ValueError("LU decomposition requires a square matrix.")

    
    L = Matrix("real", size, size)
    U = Matrix("real", size, size)

    
    L_entries = [[0.0 for _ in range(size)] for _ in range(size)]
    U_entries = [[0.0 for _ in range(size)] for _ in range(size)]

    
    for i in range(size):
        
        for j in range(i, size):
            U_entries[i][j] = matrix.entries[i][j] - sum(L_entries[i][k] * U_entries[k][j] for k in range(i))

        
        for j in range(i, size):
            if i == j:
                L_entries[i][j] = 1  
            else:
                if abs(U_entries[i][i]) < 1e-10:
                    raise ValueError("Matrix is singular and cannot be decomposed.")
                L_entries[j][i] = (matrix.entries[j][i] - sum(L_entries[j][k] * U_entries[k][i] for k in range(i))) / U_entries[i][i]

    
    L.set_entries(L_entries)
    U.set_entries(U_entries)

    return L, U

