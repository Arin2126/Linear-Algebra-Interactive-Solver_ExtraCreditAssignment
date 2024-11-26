from core_library import Matrix

def question_10():
    
    print("Welcome to Question 10: Matrix Decompositions!")
    while True:
        print("\nChoose an operation:")
        print("1. Compute Polar Decomposition of a Matrix")
        print("2. Compute Cholesky Decomposition of a Matrix")
        print("3. Compute Singular Value Decomposition (SVD) of a Matrix")
        print("4. Back to Main Menu")

        try:
            choice = int(input("Enter your choice (1-4): "))
            if choice == 4:
                print("Returning to Main Menu...")
                break

            rows = int(input("Enter the number of rows: "))
            cols = int(input("Enter the number of columns: "))
            if choice in [1, 2] and rows != cols:
                raise ValueError("Matrix must be square for Polar and Cholesky decompositions.")

            print("Enter the matrix values row by row:")
            values = [
                list(map(float, input(f"Enter row {i+1} values separated by spaces: ").split()))
                for i in range(rows)
            ]
            if any(len(row) != cols for row in values):
                raise ValueError("Each row must have the specified number of columns.")

            matrix = Matrix("real", rows, cols)
            matrix.set_entries(values)

            if choice == 1:
                U, H = compute_polar_decomposition(matrix)
                print("Unitary Matrix (U):")
                print_matrix(U)
                print("Hermitian Positive-Semidefinite Matrix (H):")
                print_matrix(H)
            elif choice == 2:
                L = compute_cholesky_decomposition(matrix)
                print("Lower Triangular Matrix (L):")
                print_matrix(L)
            elif choice == 3:
                U, Sigma, V_T = compute_svd(matrix)
                print("Left Singular Vectors (U):")
                print_matrix(U)
                print("Singular Values (Sigma):")
                print_matrix(Sigma)
                print("Right Singular Vectors Transpose (V^T):")
                print_matrix(V_T)
            else:
                print("Invalid choice. Please select a valid option.")
        except ValueError as e:
            print(f"Input Error: {e}")
        except Exception as e:
            print(f"Error: {e}")


def compute_polar_decomposition(matrix):
    """Compute the Polar Decomposition of a square matrix."""
    A_T_A = multiply_matrices(matrix.transpose(), matrix)
    H = compute_sqrt_matrix(A_T_A)  # Hermitian Positive-Semidefinite Matrix
    U = multiply_matrices(matrix, invert_matrix(H))  # Unitary Matrix
    return U, H


def compute_cholesky_decomposition(matrix):
    """Compute the Cholesky Decomposition of a Hermitian positive definite matrix."""
    size = matrix.rows
    L = Matrix("real", size, size)
    L.set_entries([[0] * size for _ in range(size)])

    for i in range(size):
        for j in range(i + 1):
            sum_k = sum(L.entries[i][k] * L.entries[j][k] for k in range(j))
            if i == j:
                L.entries[i][j] = sqrt(matrix.entries[i][i] - sum_k)
            else:
                L.entries[i][j] = (matrix.entries[i][j] - sum_k) / L.entries[j][j]
    return L


def compute_svd(matrix):
    """Compute the Singular Value Decomposition (SVD) of a matrix."""
    A_T_A = multiply_matrices(matrix.transpose(), matrix)
    eigenvalues, V = compute_eigen_decomposition(A_T_A)
    Sigma = Matrix("real", matrix.rows, matrix.cols)
    Sigma.set_entries([[sqrt(eigenvalues[i]) if i < len(eigenvalues) else 0 for i in range(matrix.cols)] for j in range(matrix.rows)])
    U = multiply_matrices(matrix, V)
    for i in range(V.cols):
        norm = sqrt(sum(U.entries[j][i] ** 2 for j in range(U.rows)))
        for j in range(U.rows):
            U.entries[j][i] /= norm
    return U, Sigma, V.transpose()

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


def compute_sqrt_matrix(matrix):
    """Compute the square root of a positive-semidefinite matrix."""
    eigenvalues, eigenvectors = compute_eigen_decomposition(matrix)
    sqrt_eigenvalues = [sqrt(val) for val in eigenvalues]
    D = Matrix("real", matrix.rows, matrix.cols)
    D.set_entries([[sqrt_eigenvalues[i] if i == j else 0 for j in range(matrix.cols)] for i in range(matrix.rows)])
    return multiply_matrices(multiply_matrices(eigenvectors, D), eigenvectors.transpose())


def compute_eigen_decomposition(matrix):
    """Compute eigenvalues and eigenvectors for a symmetric matrix."""
    size = matrix.rows
    eigenvalues = [0] * size  # Approximation of eigenvalues
    eigenvectors = Matrix("real", size, size)
    eigenvectors.set_entries([[1 if i == j else 0 for j in range(size)] for i in range(size)])
    for i in range(size):
        x = [1] * size
        for _ in range(20):  # Iterative approximation
            x = [sum(matrix.entries[row][col] * x[col] for col in range(size)) for row in range(size)]
            norm = sqrt(sum(x[i] ** 2 for i in range(size)))
            x = [xi / norm for xi in x]
        eigenvalues[i] = sum(x[row] * sum(matrix.entries[row][col] * x[col] for col in range(size)) for row in range(size))
        for j in range(size):
            eigenvectors.entries[j][i] = x[j]
    return eigenvalues, eigenvectors


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
