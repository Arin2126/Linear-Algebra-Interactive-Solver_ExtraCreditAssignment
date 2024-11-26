from core_library import Matrix, Vector, compute_rref, compute_plu_decomposition


class LinearSystem:
    
    def __init__(self, A, b):
        if A.rows != len(b.values):
            raise ValueError("Matrix A row count must match the length of vector b.")
        self.A = A
        self.b = b


def question_4():
    
    print("Welcome to Question 4: Systems of Linear Equations!")
    while True:
        print("\nChoose an operation:")
        print("1. Define a System of Linear Equations")
        print("2. Check if the System is Consistent")
        print("3. Solve the System Using Gaussian Elimination")
        print("4. Check if the Span of S1 is a Subspace of S2")
        print("5. Express the Solution Set in Terms of Free Variables")
        print("6. Solve the System Using PLU Decomposition")
        print("7. Back to Main Menu")

        try:
            choice = int(input("Enter your choice (1-7): "))
            if choice == 7:
                print("Returning to Main Menu...")
                break

            if choice == 1:
                define_linear_system()
            elif choice == 2:
                check_consistency()
            elif choice == 3:
                solve_using_gaussian_elimination()
            elif choice == 4:
                check_subspace()
            elif choice == 5:
                express_solution_with_free_variables()
            elif choice == 6:
                solve_using_plu()
            else:
                print("Invalid choice. Please select a valid option.")
        except ValueError as e:
            print(f"Input Error: {e}")
        except Exception as e:
            print(f"Error: {e}")


### Requirement (a): Define a Linear System
def define_linear_system():
    
    rows = int(input("Enter the number of rows in matrix A: "))
    cols = int(input("Enter the number of columns in matrix A: "))
    print("Enter the matrix A values row by row:")
    A_values = [
        list(map(float, input(f"Enter row {i+1} values separated by spaces: ").split()))
        for i in range(rows)
    ]
    print("Enter the elements of vector b:")
    b_values = list(map(float, input(f"Enter {rows} values separated by spaces: ").split()))

    A = Matrix("real", rows, cols)
    A.set_entries(A_values)
    b = Vector("real", rows)
    b.set_values(b_values)

    system = LinearSystem(A, b)
    print("System defined successfully!")
    return system


### Requirement (b): Check Consistency
def check_consistency():
    
    system = define_linear_system()
    augmented = system.A.append_column(system.b.values)
    rref = compute_rref(augmented)

    is_consistent = all(
        not all(abs(val) < 1e-10 for val in row[:-1]) or abs(row[-1]) < 1e-10 for row in rref.entries
    )
    if is_consistent:
        print("The system is consistent.")
    else:
        print("The system is NOT consistent.")


### Requirement (c): Solve Using Gaussian Elimination
def solve_using_gaussian_elimination():
    
    system = define_linear_system()
    augmented = system.A.append_column(system.b.values)
    rref = compute_rref(augmented)

    solutions = [row[-1] for row in rref.entries if abs(row[-1]) > 1e-10]
    print("Solution:")
    print(solutions)


### Requirement (d): Check if Span of S1 is a Subspace of S2
def check_subspace():
    
    num_vectors1 = int(input("Enter the number of vectors in S1: "))
    length1 = int(input("Enter the length of each vector in S1: "))
    S1 = []
    for i in range(num_vectors1):
        print(f"Enter the elements of vector {i+1} in S1:")
        values = list(map(float, input().split()))
        S1.append(Vector("real", length1, values))

    num_vectors2 = int(input("Enter the number of vectors in S2: "))
    length2 = int(input("Enter the length of each vector in S2: "))
    S2 = []
    for i in range(num_vectors2):
        print(f"Enter the elements of vector {i+1} in S2:")
        values = list(map(float, input().split()))
        S2.append(Vector("real", length2, values))

    if length1 != length2:
        print("The sets cannot span the same subspace as their vector lengths differ.")
        return

    combined_matrix = Matrix("real", length1, num_vectors1 + num_vectors2)
    combined_matrix.set_entries([vec.values for vec in S1 + S2])
    rref = compute_rref(combined_matrix)

    if num_vectors1 <= num_vectors2:
        print("S1 is a subspace of S2.")
    else:
        print("S1 is NOT a subspace of S2.")


### Requirement (e): Express Solution in Terms of Free Variables
def express_solution_with_free_variables():
    
    system = define_linear_system()
    augmented = system.A.append_column(system.b.values)
    rref = compute_rref(augmented)

    print("Expressing solution in terms of free variables:")
    for row in rref.entries:
        print(row)


### Requirement (f): Solve Using PLU Decomposition
def solve_using_plu():
    
    system = define_linear_system()
    P, L, U = compute_plu_decomposition(system.A)
    y = forward_substitution(L, system.b)
    x = backward_substitution(U, y)
    print("Solution using PLU decomposition:")
    print(x.values)


def forward_substitution(L, b):
    
    y = Vector("real", b.length)
    for i in range(b.length):
        y.values[i] = b.values[i] - sum(L.entries[i][j] * y.values[j] for j in range(i))
    return y


def backward_substitution(U, y):
    
    x = Vector("real", y.length)
    for i in reversed(range(y.length)):
        x.values[i] = (y.values[i] - sum(U.entries[i][j] * x.values[j] for j in range(i + 1, y.length))) / U.entries[i][i]
    return x
