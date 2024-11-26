from core_library import Matrix, Vector, compute_rref, compute_plu_decomposition, print_matrix, compute_lu_decomposition


def question_3():
    
    print("Welcome to Question 3: Elementary Operations on Matrices!")
    while True:
        print("\nChoose an operation:")
        print("1. Output Length of a Vector, Size of a Matrix, Rank, and Nullity")
        print("2. Compute RREF and Show Row Operations")
        print("3. Check if Vectors are Linearly Dependent or Independent")
        print("4. Compute Dimension of Subspace Spanned by Vectors and Find Basis")
        print("5. Compute Rank Factorization of a Matrix")
        print("6. Compute LU Decomposition of a Matrix")
        print("7. Compute PLU Decomposition of a Matrix")
        print("8. Back to Main Menu")

        try:
            choice = int(input("Enter your choice (1-8): "))
            if choice == 8:
                print("Returning to Main Menu...")
                break

            if choice == 1:
                output_matrix_properties()
            elif choice == 2:
                compute_rref_with_options()
            elif choice == 3:
                check_linear_dependence()
            elif choice == 4:
                compute_span_dimension_and_basis()
            elif choice == 5:
                compute_rank_factorization()
            elif choice == 6:
                compute_lu()
            elif choice == 7:
                compute_plu()
            else:
                print("Invalid choice. Please select a valid option.")
        except ValueError as e:
            print(f"Input Error: {e}")
        except Exception as e:
            print(f"Error: {e}")


### Requirement (a): Output Matrix Properties
def output_matrix_properties():
    
    rows = int(input("Enter the number of rows in the matrix: "))
    cols = int(input("Enter the number of columns in the matrix: "))
    print("Enter the matrix values row by row:")
    values = [
        list(map(float, input(f"Enter row {i+1} values separated by spaces: ").split()))
        for i in range(rows)
    ]
    matrix = Matrix("real", rows, cols)
    matrix.set_entries(values)

    print(f"Size of the matrix: {rows} x {cols}")
    rref = compute_rref(matrix)
    rank = sum(1 for row in rref.entries if any(abs(val) > 1e-10 for val in row))
    nullity = cols - rank
    print(f"Rank of the matrix: {rank}")
    print(f"Nullity of the matrix: {nullity}")


### Requirement (b): Compute RREF with Row Operations
def compute_rref_with_options():
    
    rows = int(input("Enter the number of rows in the matrix: "))
    cols = int(input("Enter the number of columns in the matrix: "))
    print("Enter the matrix values row by row:")
    values = [
        list(map(float, input(f"Enter row {i+1} values separated by spaces: ").split()))
        for i in range(rows)
    ]
    matrix = Matrix("real", rows, cols)
    matrix.set_entries(values)

    show_operations = input("Would you like to see row operations and elementary matrices? (yes/no): ").lower()
    if show_operations == "yes":
        print("Row operations and corresponding elementary matrices:")
        rref = compute_rref_with_steps(matrix)
    else:
        rref = compute_rref(matrix)

    print("Reduced Row Echelon Form (RREF):")
    print_matrix(rref)


def compute_rref_with_steps(matrix):
    
    rows, cols = matrix.rows, matrix.cols
    rref_matrix = [row[:] for row in matrix.entries]  

    print("Initial Matrix:")
    for row in rref_matrix:
        print(row)
    print()

    pivot_row = 0
    for pivot_col in range(cols):
        if pivot_row >= rows:
            break

        
        pivot_element_row = None
        for i in range(pivot_row, rows):
            if abs(rref_matrix[i][pivot_col]) > 1e-10:  
                pivot_element_row = i
                break

        if pivot_element_row is None:
            continue  

        
        if pivot_element_row != pivot_row:
            rref_matrix[pivot_row], rref_matrix[pivot_element_row] = rref_matrix[pivot_element_row], rref_matrix[pivot_row]
            print(f"Swapped rows {pivot_row} and {pivot_element_row}:")
            for row in rref_matrix:
                print(row)
            print()

        
        pivot_element = rref_matrix[pivot_row][pivot_col]
        for col in range(cols):
            rref_matrix[pivot_row][col] /= pivot_element
        print(f"Normalized row {pivot_row}:")
        for row in rref_matrix:
            print(row)
        print()

        
        for i in range(rows):
            if i != pivot_row:
                factor = rref_matrix[i][pivot_col]
                for col in range(cols):
                    rref_matrix[i][col] -= factor * rref_matrix[pivot_row][col]
                print(f"Eliminated column {pivot_col} in row {i}:")
                for row in rref_matrix:
                    print(row)
                print()

        
        pivot_row += 1

    print("Final RREF Matrix:")
    for row in rref_matrix:
        print(row)
    print()

    
    rref_result = Matrix(matrix.field, rows, cols)
    rref_result.set_entries(rref_matrix)
    return rref_result


### Requirement (c): Check Linear Dependence
def check_linear_dependence():
    
    num_vectors = int(input("Enter the number of vectors: "))
    length = int(input("Enter the length of each vector: "))
    vectors = []
    for i in range(num_vectors):
        print(f"Enter the elements of vector {i+1}:")
        values = list(map(float, input().split()))
        vector = Vector("real", length)
        vector.set_values(values)
        vectors.append(vector)

    matrix = Matrix("real", length, num_vectors)
    matrix.set_entries([vec.values for vec in vectors])
    rref = compute_rref(matrix)

    is_independent = all(
        any(abs(val) > 1e-10 for val in rref.entries[i]) for i in range(min(length, num_vectors))
    )
    if is_independent:
        print("The vectors are linearly independent.")
    else:
        print("The vectors are linearly dependent.")


### Requirement (d): Compute Span Dimension and Basis
def compute_span_dimension_and_basis():
    """Compute the span dimension and basis of a set of vectors."""
    try:
        # Input number of vectors and their length
        num_vectors = int(input("Enter the number of vectors: "))
        if num_vectors <= 0:
            print("The number of vectors must be a positive integer.")
            return

        vector_length = int(input("Enter the length of each vector: "))
        if vector_length <= 0:
            print("The length of each vector must be a positive integer.")
            return

        # Input the vectors
        print("Enter the vectors:")
        vectors = []
        for i in range(num_vectors):
            while True:
                try:
                    vector = list(map(float, input(f"Enter vector {i+1} values separated by spaces: ").split()))
                    if len(vector) != vector_length:
                        print(f"Error: Vector must have exactly {vector_length} values. Please try again.")
                        continue
                    vectors.append(vector)
                    break
                except ValueError:
                    print("Error: Invalid input. Please enter valid numbers separated by spaces.")

        # Perform computation
        matrix = [vector[:] for vector in vectors]
        rows, cols = len(matrix), len(matrix[0])

        pivot_cols = []
        pivot_row = 0

        for col in range(cols):
            if pivot_row >= rows:
                break

            pivot = None
            for row in range(pivot_row, rows):
                if abs(matrix[row][col]) > 1e-10:
                    pivot = row
                    break

            if pivot is None:
                continue

            if pivot != pivot_row:
                matrix[pivot], matrix[pivot_row] = matrix[pivot_row], matrix[pivot]

            pivot_element = matrix[pivot_row][col]
            for c in range(cols):
                matrix[pivot_row][c] /= pivot_element

            for row in range(rows):
                if row != pivot_row:
                    factor = matrix[row][col]
                    for c in range(cols):
                        matrix[row][c] -= factor * matrix[pivot_row][c]

            pivot_cols.append(col)
            pivot_row += 1

        basis = []
        for col in pivot_cols:
            basis_vector = [matrix[row][col] for row in range(rows)]
            basis.append(basis_vector)

        # Output results
        print(f"Dimension of the span: {len(basis)}")
        print("Basis vectors:")
        for vec in basis:
            print(vec)

    except ValueError as e:
        print(f"Input Error: {e}")
    except Exception as e:
        print(f"Unexpected Error: {e}")


### Requirement (e): Rank Factorization
def compute_rank_factorization():
    
    try:
        # Input number of rows and columns
        rows = int(input("Enter the number of rows in the matrix: "))
        cols = int(input("Enter the number of columns in the matrix: "))

        # Validate dimensions
        if rows <= 0 or cols <= 0:
            print("Matrix dimensions must be positive integers.")
            return

        # Input the matrix
        print("Enter the matrix values row by row:")
        values = []
        for i in range(rows):
            while True:
                try:
                    row = list(map(float, input(f"Enter row {i+1} values separated by spaces: ").split()))
                    if len(row) != cols:
                        print(f"Error: Row must have exactly {cols} values. Please try again.")
                        continue
                    values.append(row)
                    break
                except ValueError:
                    print("Error: Invalid input. Please enter valid numbers separated by spaces.")

        # Compute the RREF of the input matrix
        rref_matrix = [row[:] for row in values]
        pivot_rows = []
        pivot_cols = []

        pivot_row = 0
        for col in range(cols):
            if pivot_row >= rows:
                break

            pivot = None
            for row in range(pivot_row, rows):
                if abs(rref_matrix[row][col]) > 1e-10:
                    pivot = row
                    break

            if pivot is None:
                continue

            if pivot != pivot_row:
                rref_matrix[pivot], rref_matrix[pivot_row] = rref_matrix[pivot_row], rref_matrix[pivot]

            pivot_element = rref_matrix[pivot_row][col]
            for c in range(cols):
                rref_matrix[pivot_row][c] /= pivot_element

            for row in range(rows):
                if row != pivot_row:
                    factor = rref_matrix[row][col]
                    for c in range(cols):
                        rref_matrix[row][c] -= factor * rref_matrix[pivot_row][c]

            pivot_rows.append(pivot_row)
            pivot_cols.append(col)
            pivot_row += 1

        # Construct U and V matrices
        r = len(pivot_rows)
        u = [[rref_matrix[i][j] if i < r else 0 for j in pivot_cols] for i in range(rows)]
        v = [[values[pivot_rows[i]][j] if i < r else 0 for j in range(cols)] for i in range(r)]

        # Display the results
        print("Matrix U (from RREF):")
        for row in u:
            print(row)
        print("Matrix V:")
        for row in v:
            print(row)

    except ValueError as e:
        print(f"Input Error: {e}")
    except Exception as e:
        print(f"Unexpected Error: {e}")
    
    


### Requirement (f): LU Decomposition
def compute_lu():
    """Compute the LU decomposition of a matrix."""
    rows = int(input("Enter the number of rows in the matrix: "))
    cols = int(input("Enter the number of columns in the matrix: "))
    if rows != cols:
        print("LU decomposition requires a square matrix.")
        return

    print("Enter the matrix values row by row:")
    values = [
        list(map(float, input(f"Enter row {i+1} values separated by spaces: ").split()))
        for i in range(rows)
    ]
    matrix = Matrix("real", rows, cols)
    matrix.set_entries(values)

    try:
        L, U = compute_lu_decomposition(matrix)
        print("Lower Triangular Matrix (L):")
        print_matrix(L)
        print("Upper Triangular Matrix (U):")
        print_matrix(U)
    except ValueError as e:
        print(f"Error: {e}")


### Requirement (g): PLU Decomposition
def compute_plu():
    """Compute the PLU decomposition of a matrix."""
    rows = int(input("Enter the number of rows in the matrix: "))
    cols = int(input("Enter the number of columns in the matrix: "))
    if rows != cols:
        print("PLU decomposition requires a square matrix.")
        return

    print("Enter the matrix values row by row:")
    values = [
        list(map(float, input(f"Enter row {i+1} values separated by spaces: ").split()))
        for i in range(rows)
    ]
    matrix = Matrix("real", rows, cols)
    matrix.set_entries(values)

    try:
        P, L, U = compute_plu_decomposition(matrix)
        print("Permutation Matrix (P):")
        print_matrix(P)
        print("Lower Triangular Matrix (L):")
        print_matrix(L)
        print("Upper Triangular Matrix (U):")
        print_matrix(U)
    except ValueError as e:
        print(f"Error: {e}")



