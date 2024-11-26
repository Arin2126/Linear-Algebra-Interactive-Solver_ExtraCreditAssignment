from question_files.q2 import question_2
from question_files.q3 import question_3
from question_files.q4 import question_4
from question_files.q5 import question_5
from question_files.q6 import question_6
from question_files.q7 import question_7
from question_files.q8 import question_8
from question_files.q9 import question_9
from question_files.q10 import question_10

def main():
    print("Welcome to the Linear Algebra Assignment Program!")
    while True:
        
        print("\nPlease choose a question to solve:")
        print("2. Question 2: Matrix Properties")
        print("3. Question 3: Elementary Operations")
        print("4. Question 4: Solving Linear Equations")
        print("5. Question 5: Invertible Matrices")
        print("6. Question 6: Change of Basis")
        print("7. Question 7: Determinants")
        print("8. Question 8: Inner Products")
        print("9. Question 9: Eigenvalues and Eigenvectors")
        print("10. Question 10: Advanced Decompositions")
        print("11. Exit")

        try:
            
            choice = int(input("Enter your choice (2-11): "))
            
            
            if choice == 2:
                question_2()
            elif choice == 3:
                question_3()
            elif choice == 4:
                question_4()
            elif choice == 5:
                question_5()
            elif choice == 6:
                question_6()
            elif choice == 7:
                question_7()
            elif choice == 8:
                question_8()
            elif choice == 9:
                question_9()
            elif choice == 10:
                question_10()
            elif choice == 11:
                print("Thank you for using the program. Goodbye!")
                break
            else:
                print("Invalid choice. Please select a number between 2 and 11.")
        
        except ValueError:
            
            print("Invalid input. Please enter a number between 2 and 11.")
        except Exception as e:
            
            print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()
