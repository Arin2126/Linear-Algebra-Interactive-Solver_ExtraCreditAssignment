from core_library import Matrix, Vector, compute_rref, invert_matrix, sqrt, multiply_matrix_vector, multiply_matrices, print_matrix


def question_8():
    
    print("Welcome to Question 8: Inner Products!")
    while True:
        print("\nChoose an operation:")
        print("1. Compute the Inner Product of Two Vectors")
        print("2. Check if Two Vectors are Orthogonal")
        print("3. Gram-Schmidt Orthogonalization of a Set of Vectors")
        print("4. Compute QR Factorization of a Matrix")
        print("5. Compute the Moore-Penrose Pseudoinverse of a Matrix")
        print("6. Compute the Least Squares Solution of a System of Equations")
        print("7. Back to Main Menu")

        try:
            choice = int(input("Enter your choice (1-7): "))
            if choice == 7:
                print("Returning to Main Menu...")
                break

            if choice == 1:
                compute_inner_product_interactive()
            elif choice == 2:
                check_orthogonality()
            elif choice == 3:
                gram_schmidt_interactive()
            elif choice == 4:
                qr_factorization_interactive()
            elif choice == 5:
                compute_pseudoinverse_interactive()
            elif choice == 6:
                least_squares_solution_interactive()
            else:
                print("Invalid choice. Please select a valid option.")
        except ValueError as e:
            print(f"Input Error: {e}")
        except Exception as e:
            print(f"Error: {e}")




# (a) Compute the inner product of two vectors
def compute_inner_product(vector1, vector2):
    
    return sum(v1 * v2 for v1, v2 in zip(vector1.values, vector2.values))


def compute_inner_product_interactive():
    
    length = int(input("Enter the length of the vectors: "))
    print("Enter the elements of the first vector:")
    values1 = list(map(float, input().split()))
    print("Enter the elements of the second vector:")
    values2 = list(map(float, input().split()))
    vector1 = Vector("real", length)
    vector2 = Vector("real", length)
    vector1.set_values(values1)
    vector2.set_values(values2)
    print(f"Inner Product: {compute_inner_product(vector1, vector2)}")


# (b) Check if two vectors are orthogonal
def are_orthogonal(vector1, vector2):
    
    return compute_inner_product(vector1, vector2) == 0


def check_orthogonality():
    
    length = int(input("Enter the length of the vectors: "))
    print("Enter the elements of the first vector:")
    values1 = list(map(float, input().split()))
    print("Enter the elements of the second vector:")
    values2 = list(map(float, input().split()))
    vector1 = Vector("real", length)
    vector2 = Vector("real", length)
    vector1.set_values(values1)
    vector2.set_values(values2)
    print(f"Vectors are orthogonal: {are_orthogonal(vector1, vector2)}")


# (c) Gram-Schmidt orthogonalization
def gram_schmidt(vectors):
    
    orthogonalized_vectors = []
    for vector in vectors:
        projection = [sum(compute_inner_product(vector, orth_vec) * orth_vec.values[i] for orth_vec in orthogonalized_vectors) for i in range(vector.length)]
        orthogonalized = [vector.values[i] - projection[i] for i in range(vector.length)]
        norm = sqrt(sum(val ** 2 for val in orthogonalized))
        orthogonalized_vectors.append(Vector("real", vector.length, [val / norm for val in orthogonalized]))
    return orthogonalized_vectors


def gram_schmidt_interactive():
    
    num_vectors = int(input("Enter the number of vectors: "))
    length = int(input("Enter the length of each vector: "))
    vectors = []
    for i in range(num_vectors):
        print(f"Enter the elements of vector {i+1}:")
        values = list(map(float, input().split()))
        vector = Vector("real", length)
        vector.set_values(values)
        vectors.append(vector)
    orthogonalized_vectors = gram_schmidt(vectors)
    print("Orthogonalized Vectors:")
    for vec in orthogonalized_vectors:
        print(vec.values)


# (d) QR Factorization
def qr_factorization(matrix):
    
    size = matrix.rows
    Q = Matrix("real", size, size)
    R = Matrix("real", size, size)
    Q.set_entries([[0] * size for _ in range(size)])
    R.set_entries([[0] * size for _ in range(size)])
    # Gram-Schmidt process
    A_columns = [Vector("real", size, [matrix.entries[i][j] for i in range(size)]) for j in range(size)]
    Q_columns = gram_schmidt(A_columns)
    for i in range(size):
        for j in range(i, size):
            R.entries[i][j] = compute_inner_product(Q_columns[i], A_columns[j])
        for k in range(size):
            Q.entries[k][i] = Q_columns[i].values[k]
    return Q, R


def qr_factorization_interactive():
    
    size = int(input("Enter the size of the square matrix: "))
    print("Enter the matrix values row by row:")
    values = [
        list(map(float, input(f"Enter row {i+1} values separated by spaces: ").split()))
        for i in range(size)
    ]
    matrix = Matrix("real", size, size)
    matrix.set_entries(values)
    Q, R = qr_factorization(matrix)
    print("Orthogonal Matrix (Q):")
    print_matrix(Q)
    print("Upper Triangular Matrix (R):")
    print_matrix(R)


# (e) Moore-Penrose Pseudoinverse
def compute_pseudoinverse(matrix):
    
    A_T = matrix.transpose()
    A_T_A = multiply_matrices(A_T, matrix)
    A_T_A_inverse = invert_matrix(A_T_A)
    return multiply_matrices(A_T_A_inverse, A_T)


def compute_pseudoinverse_interactive():
    
    rows = int(input("Enter the number of rows: "))
    cols = int(input("Enter the number of columns: "))
    print("Enter the matrix values row by row:")
    values = [
        list(map(float, input(f"Enter row {i+1} values separated by spaces: ").split()))
        for i in range(rows)
    ]
    matrix = Matrix("real", rows, cols)
    matrix.set_entries(values)
    pseudoinverse = compute_pseudoinverse(matrix)
    print("Pseudoinverse of the Matrix:")
    print_matrix(pseudoinverse)


# (f) Least Squares Solution
def least_squares_solution(matrix, b):
    
    A_pseudoinverse = compute_pseudoinverse(matrix)
    return multiply_matrix_vector(A_pseudoinverse, b)


def least_squares_solution_interactive():
    
    rows = int(input("Enter the number of rows: "))
    cols = int(input("Enter the number of columns: "))
    print("Enter the matrix values row by row:")
    values = [
        list(map(float, input(f"Enter row {i+1} values separated by spaces: ").split()))
        for i in range(rows)
    ]
    matrix = Matrix("real", rows, cols)
    matrix.set_entries(values)
    print("Enter the right-hand side vector (b):")
    b_values = list(map(float, input(f"Enter {rows} values separated by spaces: ").split()))
    b = Vector("real", rows)
    b.set_values(b_values)
    solution = least_squares_solution(matrix, b)
    print("Least Squares Solution:")
    print(solution.values)
