from core_library import Matrix, Vector

def question_9():
    
    print("Welcome to Question 9: Eigenvalues and Eigenvectors!")
    while True:
        print("\nChoose an operation:")
        print("1. Compute Roots of a Polynomial (Aberth Method)")
        print("2. Compute Characteristic Polynomial, Minimal Polynomial, and Eigenvalues")
        print("3. Check if Two Matrices are Similar and Find Change of Basis Matrix")
        print("4. Compute Algebraic/Geometric Multiplicities and Eigen-Basis")
        print("5. Check Diagonalizability and Find Change of Basis Matrix")
        print("6. Back to Main Menu")

        try:
            choice = int(input("Enter your choice (1-6): "))
            if choice == 6:
                print("Returning to Main Menu...")
                break

            if choice == 1:
                compute_roots_of_polynomial()
            elif choice == 2:
                compute_characteristic_minimal_eigenvalues()
            elif choice == 3:
                check_matrix_similarity()
            elif choice == 4:
                compute_multiplicities_and_eigen_basis()
            elif choice == 5:
                check_diagonalizability()
            else:
                print("Invalid choice. Please select a valid option.")
        except ValueError as e:
            print(f"Input Error: {e}")
        except Exception as e:
            print(f"Error: {e}")


### Part (a): Roots of Polynomial (Aberth Method)
def compute_roots_of_polynomial():
    
    n = int(input("Enter the degree of the polynomial (n): "))
    print(f"Enter the {n + 1} coefficients of the polynomial (a_0 to a_n):")
    coefficients = list(map(complex, input().split()))
    if len(coefficients) != n + 1:
        raise ValueError("Number of coefficients must be n + 1.")

    roots = aberth_method(coefficients)
    print("Roots of the polynomial:")
    for root in roots:
        print(f"{root:.6f}")


def aberth_method(coefficients):
    
    n = len(coefficients) - 1
    roots = [complex(i, i) for i in range(1, n + 1)]  
    max_iterations = 100
    tolerance = 1e-9

    def poly(z):
        return sum(coefficients[i] * z ** i for i in range(len(coefficients)))

    def poly_derivative(z):
        return sum(i * coefficients[i] * z ** (i - 1) for i in range(1, len(coefficients)))

    for _ in range(max_iterations):
        corrections = []
        for i in range(n):
            numerator = poly(roots[i])
            denominator = poly_derivative(roots[i])
            correction = numerator / denominator
            sum_term = sum(1 / (roots[i] - roots[j]) for j in range(n) if i != j)
            correction /= (1 - correction * sum_term)
            corrections.append(correction)

        roots = [roots[i] - corrections[i] for i in range(n)]
        if max(abs(c) for c in corrections) < tolerance:
            break

    return roots


### Part (b): Characteristic Polynomial, Minimal Polynomial, and Eigenvalues
def compute_characteristic_minimal_eigenvalues():
    
    size = int(input("Enter the size of the square matrix: "))
    print("Enter the matrix values row by row:")
    values = [
        list(map(float, input(f"Enter row {i+1} values separated by spaces: ").split()))
        for i in range(size)
    ]
    matrix = Matrix("real", size, size)
    matrix.set_entries(values)

    characteristic_poly = compute_characteristic_polynomial(matrix)
    minimal_poly = compute_minimal_polynomial(matrix)
    eigenvalues = compute_eigenvalues(matrix)

    print(f"Characteristic Polynomial: {characteristic_poly}")
    print(f"Minimal Polynomial: {minimal_poly}")
    print(f"Eigenvalues: {eigenvalues}")


def compute_characteristic_polynomial(matrix):
    """Compute the characteristic polynomial of a matrix."""
    size = matrix.rows
    determinant_poly = []
    for i in range(size + 1):
        modified_matrix = Matrix("real", size, size)
        modified_matrix.set_entries([
            [matrix.entries[r][c] - (1 if r == c and i > 0 else 0) for c in range(size)]
            for r in range(size)
        ])
        determinant_poly.append(determinant(modified_matrix))
    return determinant_poly


def determinant(matrix):
    
    size = matrix.rows
    if size == 1:
        return matrix.entries[0][0]
    if size == 2:
        return matrix.entries[0][0] * matrix.entries[1][1] - matrix.entries[0][1] * matrix.entries[1][0]

    result = 0
    for col in range(size):
        sub_matrix = Matrix("real", size - 1, size - 1)
        sub_matrix.set_entries([
            [matrix.entries[i][j] for j in range(size) if j != col]
            for i in range(1, size)
        ])
        result += ((-1) ** col) * matrix.entries[0][col] * determinant(sub_matrix)
    return result


def compute_minimal_polynomial(matrix):
    
    eigenvalues = compute_eigenvalues(matrix)
    minimal_poly = [1]  # Constant term
    for eigen in eigenvalues:
        minimal_poly = multiply_polynomials(minimal_poly, [1, -eigen])
    return minimal_poly


def multiply_polynomials(poly1, poly2):
    
    degree1 = len(poly1) - 1
    degree2 = len(poly2) - 1
    result = [0] * (degree1 + degree2 + 1)
    for i in range(degree1 + 1):
        for j in range(degree2 + 1):
            result[i + j] += poly1[i] * poly2[j]
    return result


def compute_eigenvalues(matrix):
    
    characteristic_poly = compute_characteristic_polynomial(matrix)
    return aberth_method(characteristic_poly)


### Part (c): Matrix Similarity and Change of Basis Matrix
def check_matrix_similarity():
    
    size = int(input("Enter the size of the square matrices: "))
    print("Enter the first matrix values row by row:")
    values1 = [
        list(map(float, input(f"Enter row {i+1} values separated by spaces: ").split()))
        for i in range(size)
    ]
    print("Enter the second matrix values row by row:")
    values2 = [
        list(map(float, input(f"Enter row {i+1} values separated by spaces: ").split()))
        for i in range(size)
    ]

    matrix1 = Matrix("real", size, size)
    matrix1.set_entries(values1)
    matrix2 = Matrix("real", size, size)
    matrix2.set_entries(values2)

    if not are_similar_matrices(matrix1, matrix2):
        print("The matrices are not similar.")
        return

    change_of_basis = find_change_of_basis_matrix(matrix1, matrix2)
    print("The matrices are similar.")
    print("Change of Basis Matrix:")
    print_matrix(change_of_basis)


def are_similar_matrices(matrix1, matrix2):
    """Check if two matrices are similar based on their eigenvalues."""
    eigenvalues1 = compute_eigenvalues(matrix1)
    eigenvalues2 = compute_eigenvalues(matrix2)

    # Sort eigenvalues based on real and imaginary parts
    sorted_eigenvalues1 = sorted(eigenvalues1, key=lambda x: (x.real, x.imag))
    sorted_eigenvalues2 = sorted(eigenvalues2, key=lambda x: (x.real, x.imag))

    return sorted_eigenvalues1 == sorted_eigenvalues2



def find_change_of_basis_matrix(matrix1, matrix2):
    
    
    size = matrix1.rows
    identity = Matrix("real", size, size)
    identity.set_entries([[1 if i == j else 0 for j in range(size)] for i in range(size)])
    return identity
def compute_multiplicities_and_eigen_basis():
    
    size = int(input("Enter the size of the square matrix: "))
    print("Enter the matrix values row by row:")
    values = [
        list(map(float, input(f"Enter row {i+1} values separated by spaces: ").split()))
        for i in range(size)
    ]
    matrix = Matrix("real", size, size)
    matrix.set_entries(values)

    eigenvalue = float(input("Enter the eigenvalue for which to compute multiplicities and basis: "))

    algebraic_multiplicity = compute_algebraic_multiplicity(matrix, eigenvalue)
    geometric_multiplicity, eigen_basis = compute_geometric_multiplicity_and_basis(matrix, eigenvalue)

    print(f"Algebraic Multiplicity of eigenvalue {eigenvalue}: {algebraic_multiplicity}")
    print(f"Geometric Multiplicity of eigenvalue {eigenvalue}: {geometric_multiplicity}")
    print("Eigen-Basis:")
    for vec in eigen_basis:
        print(vec.values)


def compute_algebraic_multiplicity(matrix, eigenvalue):
    
    characteristic_poly = compute_characteristic_polynomial(matrix)
    return characteristic_poly.count(eigenvalue)


def compute_geometric_multiplicity_and_basis(matrix, eigenvalue):
    
    size = matrix.rows
    identity = Matrix("real", size, size)
    identity.set_entries([[1 if i == j else 0 for j in range(size)] for i in range(size)])

    # Compute (A - λI)
    lambda_matrix = Matrix("real", size, size)
    lambda_matrix.set_entries([
        [matrix.entries[i][j] - (eigenvalue if i == j else 0) for j in range(size)]
        for i in range(size)
    ])

    # Find the null space of (A - λI) to compute eigen-basis
    rref = compute_rref(lambda_matrix)
    eigen_basis = []
    for i in range(size):
        if all(val == 0 for val in rref.entries[i]):
            eigenvector = Vector("real", size)
            eigenvector.set_values([1 if i == j else 0 for j in range(size)])
            eigen_basis.append(eigenvector)

    return len(eigen_basis), eigen_basis


def compute_rref(matrix):
    """Compute the Reduced Row Echelon Form (RREF) of a matrix."""
    rref_matrix = Matrix("real", matrix.rows, matrix.cols)
    rref_matrix.set_entries([row[:] for row in matrix.entries])  # Deep copy of the matrix
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

def check_diagonalizability():
    
    size = int(input("Enter the size of the square matrix: "))
    print("Enter the matrix values row by row:")
    values = [
        list(map(float, input(f"Enter row {i+1} values separated by spaces: ").split()))
        for i in range(size)
    ]
    matrix = Matrix("real", size, size)
    matrix.set_entries(values)

    eigenvalues = compute_eigenvalues(matrix)
    total_geometric_multiplicity = 0
    eigen_basis = []

    for eigenvalue in eigenvalues:
        geometric_multiplicity, basis = compute_geometric_multiplicity_and_basis(matrix, eigenvalue)
        total_geometric_multiplicity += geometric_multiplicity
        eigen_basis.extend(basis)

    if total_geometric_multiplicity == size:
        print("The matrix is diagonalizable.")
        print("Change of Basis Matrix (columns are eigenvectors):")
        change_of_basis = Matrix("real", size, size)
        change_of_basis.set_entries([vec.values for vec in eigen_basis])
        print_matrix(change_of_basis)
    else:
        print("The matrix is not diagonalizable.")


def print_matrix(matrix):
    
    for row in matrix.entries:
        print(" ".join(f"{val:.2f}" for val in row))



