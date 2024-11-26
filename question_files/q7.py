from core_library import Matrix, compute_rref


def question_7():
    
    print("Welcome to Question 7: Determinants!")
    while True:
        print("\nChoose an operation:")
        print("1. Compute Determinant Using Cofactor Expansion")
        print("2. Compute Determinant Using PLU Decomposition")
        print("3. Compute Determinant Using RREF")
        print("4. Back to Main Menu")

        try:
            choice = int(input("Enter your choice (1-4): "))
            if choice == 4:
                print("Returning to Main Menu...")
                break

            size = int(input("Enter the size of the square matrix: "))
            print("Enter the matrix values row by row:")
            values = [
                list(map(float, input(f"Enter row {i+1} values separated by spaces: ").split()))
                for i in range(size)
            ]
            matrix = Matrix("real", size, size)
            matrix.set_entries(values)

            if choice == 1:
                determinant = cofactor_expansion(matrix)
                print(f"Determinant (Cofactor Expansion): {determinant}")
            elif choice == 2:
                determinant = plu_determinant(matrix)
                print(f"Determinant (PLU Decomposition): {determinant}")
            elif choice == 3:
                determinant = rref_determinant(matrix)
                print(f"Determinant (RREF Method): {determinant}")
            else:
                print("Invalid choice. Please select a valid option.")
        except ValueError as e:
            print(f"Input Error: {e}")
        except Exception as e:
            print(f"Error: {e}")


### Requirement (a): Cofactor Expansion
def cofactor_expansion(matrix):
    
    size = matrix.rows
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
        determinant += ((-1) ** col) * matrix.entries[0][col] * cofactor_expansion(sub_matrix)
    return determinant


### Requirement (b): PLU Decomposition
def compute_plu_decomposition(matrix):
    
    # Validate that the input is a Matrix object and square
    if not isinstance(matrix, Matrix):
        raise ValueError("Input must be a Matrix object.")
    if matrix.rows != matrix.cols:
        raise ValueError("Matrix must be square for PLU decomposition.")

    size = matrix.rows

    # Initialize permutation (P), lower triangular (L), and upper triangular (U) matrices
    P = [[1 if i == j else 0 for j in range(size)] for i in range(size)]
    L = [[0 for _ in range(size)] for _ in range(size)]
    U = [[matrix.entries[i][j] for j in range(size)] for i in range(size)]

    # Perform the decomposition
    for i in range(size):
        # Find the row with the maximum absolute value in the current column for pivoting
        max_row = i
        for k in range(i + 1, size):
            if abs(U[k][i]) > abs(U[max_row][i]):
                max_row = k

        # Swap rows in U and P for pivoting
        if max_row != i:
            U[i], U[max_row] = U[max_row], U[i]
            P[i], P[max_row] = P[max_row], P[i]

        # Check for singularity
        if U[i][i] == 0:
            raise ValueError("Matrix is singular and cannot have a PLU decomposition.")

        # Eliminate elements below the pivot in the current column
        for j in range(i + 1, size):
            factor = U[j][i] / U[i][i]
            L[j][i] = factor
            for k in range(i, size):
                U[j][k] -= factor * U[i][k]

    # Set diagonal elements of L to 1
    for i in range(size):
        L[i][i] = 1

    return P, L, U


def product_of_diagonal_elements(matrix):
    
    size = len(matrix)
    product = 1
    for i in range(size):
        product *= matrix[i][i]
    return product


def count_permutations(P):
    
    size = len(P)
    visited = [False] * size
    count = 0

    for i in range(size):
        if not visited[i]:
            j = i
            cycle_length = 0
            while not visited[j]:
                visited[j] = True
                j = P[j].index(1)  # Find the index of the 1 in row j
                cycle_length += 1
            if cycle_length > 1:
                count += cycle_length - 1

    return count


def plu_determinant(matrix):
    
    # Perform PLU decomposition
    P, L, U = compute_plu_decomposition(matrix)

    # Compute the sign of the permutation matrix
    sign = (-1) ** count_permutations(P)

    # Compute the determinant as the product of diagonal elements of U
    determinant = sign * product_of_diagonal_elements(U)
    return determinant


### Requirement (c): RREF Method
def rref_determinant(matrix):
    
    rref_matrix = compute_rref(matrix)
    scaling_factor = compute_rref_scaling_factor(matrix)
    determinant = product_of_diagonal_elements(rref_matrix) * scaling_factor
    return determinant


def compute_rref_scaling_factor(matrix):
    
    size = matrix.rows
    scaling_factor = 1
    for i in range(size):
        for j in range(size):
            if i != j and matrix.entries[i][j] != 0:
                scaling_factor *= matrix.entries[i][j]
    return scaling_factor


def prod(iterable):
    
    result = 1
    for x in iterable:
        result *= x
    return result
