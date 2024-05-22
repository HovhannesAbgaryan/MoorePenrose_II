import sympy

from src.analytical.analytical import analytical_method
from src.numerical_analytical.numerical_analytical import numerical_analytical_method

# region Main()

# region Matrix sizes

print("Enter mxn matrix sizes (m >= n)")
while True:
    while True:
        m = int(input("m = "))
        if m > 1:
            break

    while True:
        n = int(input("n = "))
        if n > 1:
            break

    if m >= n:
        break
print(f"Matrix sizes are: m = {m}, n = {n}")

# endregion Matrix sizes

# region Matrix elements

t = sympy.symbols('t', real=True)

A_t = sympy.Matrix()
print("Enter elements of A(t) matrix:")
for i in range(m):
    row = [sympy.parse_expr(element, local_dict=dict(t=t))
           for element in input().split()]
    A_t = A_t.row_insert(i, sympy.Matrix([row]))

# endregion Matrix elements

# region A1(t), A2(t)

A1_t, A2_t = A_t.as_real_imag()
print(f"A1(t) = {A1_t}")
print(f"A2(t) = {A2_t}")

# endregion A1(t), A2(t)

# region Analytical Method

print("ANALYTICAL METHOD")
analytical_method(A1_t, A2_t)

# endregion Analytical Method

# region Numerical-analytical method

print("NUMERICAL-ANALYTICAL METHOD")
numerical_analytical_method(A1_t, A2_t, t)

# endregion Numerical-analytical method

# endregion Main()
