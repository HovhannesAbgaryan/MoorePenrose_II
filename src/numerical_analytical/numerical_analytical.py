import numpy
import sympy
import math


# region Functions

def find_derivatives(matrix, symbol):
    derivatives = list([matrix])
    while not matrix.is_zero_matrix:
        derivative = sympy.diff(matrix, symbol)
        derivatives.append(derivative)
        matrix = derivative
    return derivatives


def compute_discretes(derivatives, symbol, scale_coeff, value):
    discretes = list()
    for k in range(len(derivatives)):
        discrete = (scale_coeff ** k / math.factorial(k) *
                    derivatives[k].subs(symbol, value))
        discretes.append(discrete)
    return discretes


# endregion Functions

def numerical_analytical_method(real_matrix, imaginary_matrix, variable):
    m = sympy.shape(real_matrix)[0]

    # region Derivatives and discretes

    real_matrix_derives = find_derivatives(real_matrix, variable)
    imaginary_matrix_derives = find_derivatives(imaginary_matrix, variable)

    H = int(input("Enter scaling coefficient: H = "))
    appr_cntr = int(input("Enter approximation center's value: "))
    real_matrix_discretes = compute_discretes(real_matrix_derives, variable, H, appr_cntr)
    imaginary_matrix_discretes = compute_discretes(imaginary_matrix_derives, variable, H, appr_cntr)

    # endregion Derivatives and discretes

    # region C1(K), C2(K), C(K)

    C1_K = list()
    C2_K = list()
    K = max(len(real_matrix_discretes), len(imaginary_matrix_discretes))
    for k in range(K + 1):
        C1_k = numpy.zeros(shape=(m, m), dtype=numpy.int_)
        C2_k = numpy.zeros(shape=(m, m), dtype=numpy.int_)
        for i in range(k + 1):
            index1_1 = i if i < len(real_matrix_discretes) else len(real_matrix_discretes) - 1
            index1_2 = k - i if k - i < len(real_matrix_discretes) else len(real_matrix_discretes) - 1
            index2_1 = i if i < len(imaginary_matrix_discretes) else len(imaginary_matrix_discretes) - 1
            index2_2 = k - i if k - i < len(imaginary_matrix_discretes) else len(imaginary_matrix_discretes) - 1
            C1_k += (real_matrix_discretes[index1_1] * real_matrix_discretes[index1_2].T +
                     imaginary_matrix_discretes[index2_1] * imaginary_matrix_discretes[index2_2].T)
            C2_k += (imaginary_matrix_discretes[index2_1] * real_matrix_discretes[index1_2].T -
                     real_matrix_discretes[index1_1] * imaginary_matrix_discretes[index2_2].T)
        C1_K.append(C1_k)
        C2_K.append(C2_k)
    print(f"C1(K) = {C1_K}")
    print(f"C2(K) = {C2_K}")

    C_K = list(())
    for k in range(K + 1):
        C_K.append([
            [C1_K[k], -C2_K[k]],
            [C2_K[k], C1_K[k]]
        ])

    # endregion C1(K), C2(K), C(K)

    # region Y(K), Y1(K), Y2(K)

    Y_K = list(([
                    [C1_K[0].pinv(), -C2_K[0].pinv()],
                    [C2_K[0].pinv(), C1_K[0].pinv()]
                ],
    ))
    for k in range(1, K + 1):
        summa = [
            [numpy.zeros(shape=(m, m), dtype=numpy.int_), numpy.zeros(shape=(m, m), dtype=numpy.int_)],
            [numpy.zeros(shape=(m, m), dtype=numpy.int_), numpy.zeros(shape=(m, m), dtype=numpy.int_)]
        ]
        for i in range(1, k + 1):
            summa[0][0] += C_K[i][0][0] * Y_K[k - i][0][0] + C_K[i][0][1] * Y_K[k - i][1][0]
            summa[0][1] += C_K[i][0][0] * Y_K[k - i][0][1] + C_K[i][0][1] * Y_K[k - i][1][1]
            summa[1][0] += C_K[i][1][0] * Y_K[k - i][0][0] + C_K[i][1][1] * Y_K[k - i][1][0]
            summa[1][1] += C_K[i][1][0] * Y_K[k - i][0][1] + C_K[i][1][1] * Y_K[k - i][1][1]
        Y_K.append(
            [
                [-C1_K[0].pinv() * summa[0][0] + C2_K[0].pinv() * summa[1][0],
                 -C1_K[0].pinv() * summa[0][1] + C2_K[0].pinv() * summa[1][1]],
                [-C2_K[0].pinv() * summa[0][0] + -C1_K[0].pinv() * summa[1][0],
                 -C2_K[0].pinv() * summa[0][1] + -C1_K[0].pinv() * summa[1][1]]
            ])
    # print(f"Y(K) = {Y_K}")

    Y1_K = list()
    Y2_K = list()
    for k in range(K + 1):
        Y1_K.append(Y_K[k][0][0])
        Y2_K.append(Y_K[k][1][0])
    print(f"Y1(K) = {Y1_K}")
    print(f"Y2(K) = {Y2_K}")

    # endregion Y(K), Y1(K), Y2(K)

    # region Y1(t), Y2(t)

    Y1_t = numpy.zeros(shape=(m, m), dtype=numpy.int_)
    Y2_t = numpy.zeros(shape=(m, m), dtype=numpy.int_)
    for k in range(K + 1):
        Y1_t += ((variable - appr_cntr) / H) ** k * Y1_K[k]
        Y2_t += ((variable - appr_cntr) / H) ** k * Y2_K[k]
    print(f"Y1(t) = {Y1_t}")
    print(f"Y2(t) = {Y2_t}")

    # endregion Y1(t), Y2(t)

    # region A+(t)

    a = real_matrix - sympy.I * imaginary_matrix
    y = Y1_t + sympy.I * Y2_t
    A_pseudo_inverse = sympy.expand(a.T * y)
    print(f"A+(t) = {A_pseudo_inverse}")

    # endregion A+(t)
