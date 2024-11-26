from core_library import Matrix, compute_rref, compute_lu_decomposition, cofactor_expansion, multiply_matrices


def question_2():
    
    print("Welcome to Question 2: Matrix Properties!")
    while True:
        print("\nChoose a property to check:")
        print("1. Zero Matrix")
        print("2. Symmetric Matrix")
        print("3. Hermitian Matrix")
        print("4. Square Matrix")
        print("5. Orthogonal Matrix")
        print("6. Unitary Matrix")
        print("7. Scalar Matrix")
        print("8. Singular Matrix")
        print("9. Invertible Matrix")
        print("10. Identity Matrix")
        print("11. Nilpotent Matrix")
        print("12. Diagonalizable Matrix")
        print("13. LU Decomposition")
        print("14. Back to Main Menu")

        try:
            choice = int(input("Enter your choice (1-14): "))
            if choice == 14:
                print("Returning to Main Menu...")
                break

            rows = int(input("Enter the number of rows: "))
            cols = int(input("Enter the number of columns: "))
            print("Enter the matrix values row by row:")
            values = [
                list(map(float, input(f"Enter row {i+1} values separated by spaces: ").split()))
                for i in range(rows)
            ]
            matrix = Matrix("real", rows, cols)
            matrix.set_entries(values)

            if choice == 1:
                print(f"Is zero matrix: {is_zero_matrix(matrix)}")
            elif choice == 2:
                print(f"Is symmetric matrix: {is_symmetric(matrix)}")
            elif choice == 3:
                print(f"Is hermitian matrix: {is_hermitian(matrix)}")
            elif choice == 4:
                print(f"Is square matrix: {is_square(matrix)}")
            elif choice == 5:
                print(f"Is orthogonal matrix: {is_orthogonal(matrix)}")
            elif choice == 6:
                print(f"Is unitary matrix: {is_unitary(matrix)}")
            elif choice == 7:
                print(f"Is scalar matrix: {is_scalar(matrix)}")
            elif choice == 8:
                print(f"Is singular matrix: {is_singular(matrix)}")
            elif choice == 9:
                print(f"Is invertible matrix: {is_invertible(matrix)}")
            elif choice == 10:
                print(f"Is identity matrix: {is_identity(matrix)}")
            elif choice == 11:
                print(f"Is nilpotent matrix: {is_nilpotent(matrix)}")
            elif choice == 12:
                print(f"Is diagonalizable matrix: {is_diagonalizable(matrix)}")
            elif choice == 13:
                check_lu_decomposition(matrix)
            else:
                print("Invalid choice. Please select a valid option.")
        except ValueError as e:
            print(f"Input Error: {e}")
        except Exception as e:
            print(f"Error: {e}")




def is_zero_matrix(matrix):
    """Check if the matrix is a zero matrix."""
    return all(all(val == 0 for val in row) for row in matrix.entries)


def is_symmetric(matrix):
    """Check if the matrix is symmetric."""
    if not is_square(matrix):
        return False
    return all(matrix.entries[i][j] == matrix.entries[j][i] for i in range(matrix.rows) for j in range(matrix.cols))


def is_hermitian(matrix):
    """Check if the matrix is hermitian."""
    if not is_square(matrix):
        return False
    return all(matrix.entries[i][j] == matrix.entries[j][i] for i in range(matrix.rows) for j in range(matrix.cols))


def is_square(matrix):
    """Check if the matrix is square."""
    return matrix.rows == matrix.cols


def is_orthogonal(matrix):
    """Check if the matrix is orthogonal."""
    if not is_square(matrix):
        return False
    transpose = matrix.transpose()
    product = multiply_matrices(matrix, transpose)
    return is_identity(product)


def is_unitary(matrix):
    """Check if the matrix is unitary."""
    if not is_square(matrix):
        return False
    conjugate_transpose = matrix.transpose()  
    product = multiply_matrices(matrix, conjugate_transpose)
    return is_identity(product)


def is_scalar(matrix):
    """Check if the matrix is scalar."""
    if not is_square(matrix):
        return False
    diagonal_value = matrix.entries[0][0]
    return all(matrix.entries[i][j] == (diagonal_value if i == j else 0) for i in range(matrix.rows) for j in range(matrix.cols))


def is_singular(matrix):
    """Check if the matrix is singular."""
    determinant = cofactor_expansion(matrix)
    return abs(determinant) < 1e-10


def is_invertible(matrix):
    """Check if the matrix is invertible."""
    return not is_singular(matrix)


def is_identity(matrix):
    """Check if the matrix is an identity matrix."""
    if not is_square(matrix):
        return False
    return all(matrix.entries[i][j] == (1 if i == j else 0) for i in range(matrix.rows) for j in range(matrix.cols))


def is_nilpotent(matrix):
    """Check if the matrix is nilpotent."""
    if not is_square(matrix):
        return False
    result = matrix
    for _ in range(matrix.rows):
        result = multiply_matrices(result, matrix)
        if is_zero_matrix(result):
            return True
    return False


def is_diagonalizable(matrix):
    """Check if the matrix is diagonalizable."""
    if not is_square(matrix):
        return False
    rref = compute_rref(matrix)
    eigen_count = sum(1 for row in rref.entries if any(abs(val) > 1e-10 for val in row))
    return eigen_count == matrix.rows


def check_lu_decomposition(matrix):
    """Check if the matrix has LU decomposition."""
    try:
        compute_lu_decomposition(matrix)
        print("LU decomposition is possible.")
    except ValueError:
        print("LU decomposition is NOT possible.")
