import sympy


def analytical_method(real_matrix, imaginary_matrix):
    # region C1(t), C2(t), C(t)

    C1_t = real_matrix.T * real_matrix + imaginary_matrix.T * imaginary_matrix
    C2_t = real_matrix.T * imaginary_matrix - imaginary_matrix.T * real_matrix
    C_t = C1_t + C2_t * sympy.I

    # endregion C1(t), C2(t), C(t)

    # region Y(t), Y1(t), Y2(t)

    Y_t = sympy.expand(C_t.pinv())
    Y1_t, Y2_t = Y_t.as_real_imag()

    # endregion Y(t), Y1(t), Y2(t)

    # region A+(t)

    real = Y1_t * real_matrix.T + Y2_t * imaginary_matrix.T
    imaginary = sympy.expand(Y2_t * real_matrix.T - Y1_t * imaginary_matrix.T)
    A_pseudo_inverse = real + imaginary * sympy.I
    print(f"A+(t) = {A_pseudo_inverse}")

    # endregion A+(t)
