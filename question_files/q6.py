from core_library import Matrix, Vector, compute_rref, invert_matrix, multiply_matrices, print_matrix



def question_6():
    
    print("Welcome to Question 6: Coordinates and Change of Basis!")
    while True:
        print("\nChoose an operation:")
        print("1. Check if a Vector is in the Linear Span of a Set of Vectors")
        print("2. Represent a Vector as a Linear Combination of Vectors in a Set")
        print("3. Check if Two Sets Span the Same Subspace")
        print("4. Compute Coordinates of a Vector in an Ordered Basis")
        print("5. Compute the Change of Basis Matrix")
        print("6. Convert Coordinates Between Two Ordered Bases")
        print("7. Back to Main Menu")

        try:
            choice = int(input("Enter your choice (1-7): "))
            if choice == 7:
                print("Returning to Main Menu...")
                break

            if choice == 1:
                check_vector_in_span()
            elif choice == 2:
                represent_vector_as_linear_combination()
            elif choice == 3:
                check_same_subspace()
            elif choice == 4:
                compute_coordinates_in_basis()
            elif choice == 5:
                compute_change_of_basis_matrix()
            elif choice == 6:
                convert_coordinates_between_bases()
            else:
                print("Invalid choice. Please select a valid option.")
        except ValueError as e:
            print(f"Input Error: {e}")
        except Exception as e:
            print(f"Error: {e}")


### Requirement (a): Check if a Vector is in the Linear Span of a Set of Vectors
def check_vector_in_span():
    """Check if a vector is in the span of a set of vectors."""
    num_vectors = int(input("Enter the number of vectors in the set: "))
    length = int(input("Enter the length of each vector: "))
    vectors = []
    for i in range(num_vectors):
        print(f"Enter the elements of vector {i+1}:")
        values = list(map(float, input().split()))
        vector = Vector("real", length)
        vector.set_values(values)
        vectors.append(vector)
    print("Enter the elements of the vector to check:")
    v_values = list(map(float, input().split()))
    v = Vector("real", length)
    v.set_values(v_values)

    augmented_matrix = Matrix("real", length, num_vectors + 1)
    augmented_matrix.set_entries([vec.values + [v.values[i]] for i, vec in enumerate(vectors)])

    rref = compute_rref(augmented_matrix)
    if all(abs(rref.entries[i][-1]) < 1e-10 for i in range(length)):
        print("The vector is in the span of the set.")
    else:
        print("The vector is NOT in the span of the set.")


### Requirement (b): Represent a Vector as a Linear Combination of Vectors in a Set
def represent_vector_as_linear_combination():
    """Represent a vector as a linear combination of a set of vectors."""
    num_vectors = int(input("Enter the number of vectors in the set: "))
    length = int(input("Enter the length of each vector: "))
    vectors = []
    for i in range(num_vectors):
        print(f"Enter the elements of vector {i+1}:")
        values = list(map(float, input().split()))
        vector = Vector("real", length)
        vector.set_values(values)
        vectors.append(vector)
    print("Enter the elements of the vector to represent:")
    v_values = list(map(float, input().split()))
    v = Vector("real", length)
    v.set_values(v_values)

    matrix = Matrix("real", length, num_vectors)
    matrix.set_entries([vec.values for vec in vectors])
    solution = compute_rref(matrix.append_column(v.values))
    coefficients = [solution.entries[i][-1] for i in range(num_vectors)]

    print("The vector can be represented as:")
    print(f" + ".join(f"{coefficients[i]} * Vector{i+1}" for i in range(num_vectors)))


### Requirement (c): Check if Two Sets Span the Same Subspace
def check_same_subspace():
    """Check if two sets of vectors span the same subspace."""
    num_vectors1 = int(input("Enter the number of vectors in the first set: "))
    length1 = int(input("Enter the length of each vector in the first set: "))
    set1 = []
    for i in range(num_vectors1):
        print(f"Enter the elements of vector {i+1} in the first set:")
        values = list(map(float, input().split()))
        vector = Vector("real", length1)
        vector.set_values(values)
        set1.append(vector)

    num_vectors2 = int(input("Enter the number of vectors in the second set: "))
    length2 = int(input("Enter the length of each vector in the second set: "))
    set2 = []
    for i in range(num_vectors2):
        print(f"Enter the elements of vector {i+1} in the second set:")
        values = list(map(float, input().split()))
        vector = Vector("real", length2)
        vector.set_values(values)
        set2.append(vector)

    if length1 != length2:
        print("The sets cannot span the same subspace as their vector lengths differ.")
        return

    combined_set = Matrix("real", length1, num_vectors1 + num_vectors2)
    combined_set.set_entries([vec.values for vec in set1 + set2])

    rref = compute_rref(combined_set)
    rank1 = sum(1 for row in rref.entries if any(abs(val) > 1e-10 for val in row[:num_vectors1]))
    rank2 = sum(1 for row in rref.entries if any(abs(val) > 1e-10 for val in row[num_vectors1:]))

    if rank1 == rank2 == length1:
        print("The two sets span the same subspace.")
    else:
        print("The two sets do NOT span the same subspace.")

### Requirement (d): Compute Coordinates in Basis and Reconstruct Vector
def compute_coordinates_in_basis():
    """Compute coordinates of a vector in a basis and reconstruct from coordinates."""
    num_vectors = int(input("Enter the number of basis vectors: "))
    length = int(input("Enter the length of each vector: "))
    basis = []
    for i in range(num_vectors):
        print(f"Enter the elements of basis vector {i+1}:")
        values = list(map(float, input().split()))
        vector = Vector("real", length)
        vector.set_values(values)
        basis.append(vector)
    print("Enter the elements of the vector to compute coordinates:")
    v_values = list(map(float, input().split()))
    v = Vector("real", length)
    v.set_values(v_values)

    basis_matrix = Matrix("real", length, num_vectors)
    basis_matrix.set_entries([vec.values for vec in basis])
    solution = compute_rref(basis_matrix.append_column(v.values))
    coordinates = [solution.entries[i][-1] for i in range(num_vectors)]

    print(f"Coordinates of the vector in the basis: {coordinates}")

    print("Reconstructed Vector:")
    reconstructed = Vector("real", length)
    reconstructed.set_values([sum(coordinates[j] * basis[j].values[i] for j in range(num_vectors)) for i in range(length)])
    print(reconstructed.values)


### Requirement (e): Compute Change of Basis Matrix
def compute_change_of_basis_matrix():
    """Compute the change of basis matrix between two bases."""
    num_vectors = int(input("Enter the number of basis vectors in each basis: "))
    length = int(input("Enter the length of each vector: "))
    basis1 = []
    basis2 = []

    for i in range(num_vectors):
        print(f"Enter the elements of basis vector {i+1} in the first basis:")
        values = list(map(float, input().split()))
        vector = Vector("real", length)
        vector.set_values(values)
        basis1.append(vector)

    for i in range(num_vectors):
        print(f"Enter the elements of basis vector {i+1} in the second basis:")
        values = list(map(float, input().split()))
        vector = Vector("real", length)
        vector.set_values(values)
        basis2.append(vector)

    basis_matrix1 = Matrix("real", length, num_vectors)
    basis_matrix1.set_entries([vec.values for vec in basis1])

    basis_matrix2 = Matrix("real", length, num_vectors)
    basis_matrix2.set_entries([vec.values for vec in basis2])

    change_of_basis_matrix = multiply_matrices(invert_matrix(basis_matrix1), basis_matrix2)
    print("Change of Basis Matrix:")
    print_matrix(change_of_basis_matrix)


### Requirement (f): Convert Coordinates Between Two Bases
def convert_coordinates_between_bases(vector, basis1, basis2):
    n = len(basis1)
    m = len(basis2)

    if n != len(vector) or any(len(b) != n for b in basis1) or any(len(b) != m for b in basis2):
        raise ValueError("Incompatible dimensions for vector or bases.")

    matrix_basis1 = [list(b) for b in basis1]
    matrix_basis2 = [list(b) for b in basis2]

    inverse_basis1 = invert_matrix(matrix_basis1)
    if inverse_basis1 is None:
        raise ValueError("Basis1 is not invertible.")

    coordinates_in_basis1 = [sum(inverse_basis1[i][j] * vector[j] for j in range(n)) for i in range(n)]

    coordinates_in_basis2 = [sum(matrix_basis2[j][i] * coordinates_in_basis1[j] for j in range(m)) for i in range(m)]

    return coordinates_in_basis2

   



