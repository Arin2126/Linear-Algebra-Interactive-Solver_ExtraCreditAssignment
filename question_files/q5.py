from core_library import Matrix, invert_matrix, cofactor_expansion


def question_5():
    
    print("Welcome to Question 5: Invertible Matrices!")
    while True:
        print("\nChoose an operation:")
        print("1. Check if a Matrix is Square and Invertible")
        print("2. Compute the Inverse of a Matrix by Row Reduction")
        print("3. Compute the Inverse of a Matrix by Adjoint Method")
        print("4. Back to Main Menu")

        try:
            choice = int(input("Enter your choice (1-4): "))
            if choice == 4:
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
                check_square_and_invertible(matrix)
            elif choice == 2:
                compute_inverse_by_row_reduction(matrix)
            elif choice == 3:
                compute_inverse_by_adjoint(matrix)
            else:
                print("Invalid choice. Please select a valid option.")
        except ValueError as e:
            print(f"Input Error: {e}")
        except Exception as e:
            print(f"Error: {e}")


### Requirement (a): Check if a Matrix is Square and Invertible
def check_square_and_invertible(matrix):
    
    if matrix.rows != matrix.cols:
        print("The matrix is not square and therefore not invertible.")
        return

    determinant = cofactor_expansion(matrix)
    if abs(determinant) < 1e-10:
        print("The matrix is square but NOT invertible (determinant is 0).")
    else:
        print("The matrix is square and invertible (non-zero determinant).")


### Requirement (b): Compute the Inverse Using Row Reduction
def compute_inverse_by_row_reduction(matrix):
    
    if matrix.rows != matrix.cols:
        print("The matrix is not square and therefore not invertible.")
        return

    determinant = cofactor_expansion(matrix)
    if abs(determinant) < 1e-10:
        print("The matrix is NOT invertible (determinant is 0).")
        return

    try:
        inverse_matrix = invert_matrix(matrix)
        print("Inverse Matrix (Row Reduction):")
        print_matrix(inverse_matrix)
    except Exception as e:
        print(f"Error: {e}")


### Requirement (c): Compute the Inverse Using Adjoint Method
def compute_inverse_by_adjoint(matrix):
    
    if matrix.rows != matrix.cols:
        print("The matrix is not square and therefore not invertible.")
        return

    determinant = cofactor_expansion(matrix)
    if abs(determinant) < 1e-10:
        print("The matrix is NOT invertible (determinant is 0).")
        return

    adjoint_matrix = compute_adjoint(matrix)
    inverse_matrix = multiply_matrix_by_scalar(adjoint_matrix, 1 / determinant)
    print("Inverse Matrix (Adjoint Method):")
    print_matrix(inverse_matrix)


def compute_adjoint(matrix):
    
    size = matrix.rows
    adjoint_matrix = Matrix("real", size, size)
    adjoint_entries = [[0] * size for _ in range(size)]

    for i in range(size):
        for j in range(size):
            sub_matrix = Matrix("real", size - 1, size - 1)
            sub_matrix.set_entries([
                [matrix.entries[x][y] for y in range(size) if y != j]
                for x in range(size) if x != i
            ])
            cofactor = ((-1) ** (i + j)) * cofactor_expansion(sub_matrix)
            adjoint_entries[j][i] = cofactor  

    adjoint_matrix.set_entries(adjoint_entries)
    return adjoint_matrix


def multiply_matrix_by_scalar(matrix, scalar):
    
    scaled_matrix = Matrix("real", matrix.rows, matrix.cols)
    scaled_matrix.set_entries([[val * scalar for val in row] for row in matrix.entries])
    return scaled_matrix


def print_matrix(matrix):
    
    for row in matrix.entries:
        print(" ".join(f"{val:.2f}" for val in row))
