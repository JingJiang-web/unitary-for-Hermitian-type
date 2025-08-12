import numpy as np
from fractions import Fraction

a1 = np.array([0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, 0.5])
a2 = np.array([1, 1, 0, 0, 0, 0, 0, 0])
a3 = np.array([-1, 1, 0, 0, 0, 0, 0, 0])
a4 = np.array([0, -1, 1, 0, 0, 0, 0, 0])
a5 = np.array([0, 0, -1, 1, 0, 0, 0, 0])
a6 = np.array([0, 0, 0, -1, 1, 0, 0, 0])
a7 = np.array([0, 0, 0, 0, -1, 1, 0, 0])

rho_e6 = np.array([0, 1, 2, 3, 4, -4, -4, 4])
rho_e7 = np.array([0, 1, 2, 3, 4, 5, -8.5, 8.5])

beta_e6 = np.array([0.5, 0.5, 0.5, 0.5, 0.5, -0.5, -0.5, 0.5])
beta_e7 = np.array([0, 0, 0, 0, 0, 0, -1, 1])

xi_e6 = np.array([0, 0, 0, 0, 0, -2 / 3, -2 / 3, 2 / 3])
xi_e7 = np.array([0, 0, 0, 0, 0, 1, -0.5, 0.5])


def reflection(a, x):
    return x - 2 * np.dot(a, x) / np.dot(a, a) * a
def s1(x): return reflection(a1, x)
def s2(x): return reflection(a2, x)
def s3(x): return reflection(a3, x)
def s4(x): return reflection(a4, x)
def s5(x): return reflection(a5, x)
def s6(x): return reflection(a6, x)
def s7(x): return reflection(a7, x)

def compose(*functions):
    def composer(value):
        for function in reversed(functions):
            value = function(value)
        return value

    return composer


function_dict = {i: globals()[f's{i}'] for i in range(1, 8)}


def process_input(n, k):
    if n == 6:
        coxeter_type = "E6"
    elif n == 7:
        coxeter_type = "E7"
    else:
        raise ValueError("Unsupported value of n. Only 6 and 7 are supported.")
    return coxeter_type


n = int(input("n = "))
k = int(input("k = "))
w_input = input("w (e.g., [1, 3, 4]): ")
# 处理输入的 w 字符串为列表
w = [int(i) for i in w_input.strip('[]').split(',') if i.strip()]

coxeter_type = process_input(n, k)
results = []
filtered_results = []

if coxeter_type == "E6":
    composed_function = compose(*(function_dict[i] for i in w))
    result = -1 * composed_function(rho_e6)
    inner_product = np.dot(result, beta_e6)
    results.append(inner_product)

    if k == 0:
        if (inner_product == 11):
            print(f"-wrho-rho is {result - rho_e6}")
            print(f"L_w is unitary")
        else:
            print(f"-wrho-rho is {result-rho_e6}")
            print(f"L_w is not unitary")


    elif k == 1 and inner_product == 8:
        filtered_results.append((1, result, inner_product))
        lambda_0_e6 = result - rho_e6 - inner_product * xi_e6
        print(f"-wrho-rho = {(result-rho_e6).tolist()}")
        #print(f"lambda_0 is {lambda_0_e6}")
        print(f"L_w is unitary")
    elif k == 2:
        lambda_0_e6 = result - rho_e6 - inner_product * xi_e6
        lambda_0_e6_fractions = [Fraction(str(num)).limit_denominator() for num in lambda_0_e6]
        str_arr = []
        for fraction in lambda_0_e6_fractions:
            numerator = fraction.numerator
            denominator = fraction.denominator
            if denominator == 1:
                str_fraction = str(numerator)
            elif numerator < 0 and denominator < 0:
                str_fraction = f"{abs(numerator)}/{abs(denominator)}"
            elif numerator < 0:
                str_fraction = f"-{abs(numerator)}/{denominator}"
            else:
                str_fraction = f"{numerator}/{denominator}"
            str_arr.append(str_fraction)

        conditiona1 = (
                -1 * lambda_0_e6_fractions[0] == lambda_0_e6_fractions[1] == lambda_0_e6_fractions[2] ==
                lambda_0_e6_fractions[3] == lambda_0_e6_fractions[4]
                and (inner_product <= 5) and np.allclose(lambda_0_e6_fractions[0] % 1, 0.5)
        )

        conditiona2 = (
                -1 * lambda_0_e6_fractions[0] == lambda_0_e6_fractions[1] == lambda_0_e6_fractions[2] ==
                lambda_0_e6_fractions[3]
                and (inner_product <= 4)
        )

        conditiona3 = (
                lambda_0_e6_fractions[0] == lambda_0_e6_fractions[1] == lambda_0_e6_fractions[2] == 0
                and (inner_product <= 4)
        )

        conditiona4 = (
                -1 * lambda_0_e6_fractions[0] == lambda_0_e6_fractions[1] == lambda_0_e6_fractions[2]
                and (inner_product <= 3)
        )

        conditiona5 = (
                -1 * lambda_0_e6_fractions[0] == lambda_0_e6_fractions[1]
                and (inner_product <= 2)
        )

        conditiona6 = (lambda_0_e6[0] + lambda_0_e6[1] > 0
                       and (inner_product <= 1)
                       )

        condition2 = (
                lambda_0_e6_fractions[0] == lambda_0_e6_fractions[1] == lambda_0_e6_fractions[2] ==
                lambda_0_e6_fractions[3] == 0 < lambda_0_e6_fractions[4]
                and (inner_product <= 4 or inner_product == 7) and float(lambda_0_e6_fractions[4]) - 3 * float(
            lambda_0_e6_fractions[5]) == -22
        )

        condition3 = (
                lambda_0_e6_fractions[0] == lambda_0_e6_fractions[1] == lambda_0_e6_fractions[2] ==
                lambda_0_e6_fractions[3] == lambda_0_e6_fractions[4] == 0
                and (inner_product <= 8 or inner_product == 11) and float(lambda_0_e6_fractions[5]) == 22 / 3
        )

        if (
                conditiona1 or conditiona2 or conditiona3 or conditiona4 or conditiona5 or conditiona6 or condition2 or condition3):
            filtered_results.append((1, result, inner_product))
            #print(f"-wrho is {result}")
            print(f"-wrho-rho = {(result-rho_e6).tolist()}")
            #print(f"lambda_0 is {lambda_0_e6}")
            print(f"L_w is unitary")
        else:
            #print(f"-wrho is {result}")
            print(f"-wrho-rho = {(result-rho_e6).tolist()}")
            #print(f"lambda_0 is {lambda_0_e6}")
            print(f"L_w is not unitary")
    elif k >= 3:
        print(f"Error")
    else:
        #print(f"-wrho is {result}")
        lambda_0_e6 = result - rho_e6 - inner_product * xi_e6
        print(f"-wrho-rho = {(result-rho_e6).tolist()}")
        #print(f"lambda_0 is {lambda_0_e6}")
        print(f"L_w is not unitary")


elif coxeter_type == "E7":
    composed_function = compose(*(function_dict[i] for i in w))
    result = -1 * composed_function(rho_e7)
    inner_product = np.dot(result, beta_e7)
    results.append(inner_product)

    if k == 0:
        if (inner_product == 17):
            print(f"-wrho-rho = {(result-rho_e7).tolist()}")
            print(f"L_w is unitary")
        else:
            print(f"-wrho-rho = {(result-rho_e7).tolist()}")
            print(f"L_w is not unitary")

    elif (k == 1 and inner_product == 13) or (k == 2 and inner_product == 9):
        lambda_0_e7 = result - rho_e7 - inner_product * xi_e7
        condition1 = (
                0 < lambda_0_e7[0] == lambda_0_e7[1] == lambda_0_e7[2] == lambda_0_e7[3] == lambda_0_e7[4]
                and np.allclose(lambda_0_e7[0] % 1, 0.5))
        condition2 = (
                lambda_0_e7[0] == lambda_0_e7[1] == lambda_0_e7[2] == lambda_0_e7[3] == 0 < lambda_0_e7[4]
                and np.allclose(lambda_0_e7[0] % 1, 0))
        condition3 = (
                lambda_0_e7[0] == lambda_0_e7[1] == lambda_0_e7[2] == lambda_0_e7[3] == lambda_0_e7[4] == 0)
        if condition1 or condition2 or condition3:
            filtered_results.append((1, result, inner_product))
            lambda_0_e7 = result - rho_e7 - inner_product * xi_e7
            print(f"-wrho-rho = {(result-rho_e7).tolist()}")
            #print(f"lambda_0 is {lambda_0_e7}")
            print(f"L_w is unitary")
        else:
            print(f"-wrho-rho = {(result-rho_e7).tolist()}")
            #print(f"lambda_0 is {lambda_0_e7}")
            print(f"L_w is not unitary")
    elif k == 3:
        lambda_0_e7 = result - rho_e7 - inner_product * xi_e7
        conditiona1 = (
                0 < lambda_0_e7[0] == lambda_0_e7[1] == lambda_0_e7[2] == lambda_0_e7[3] == lambda_0_e7[4]
                and (inner_product <= 6) and np.allclose(lambda_0_e7[0] % 1, 0.5)
        )

        conditiona2 = (
                lambda_0_e7[0] == lambda_0_e7[1] == lambda_0_e7[2] == 0 < lambda_0_e7[3]
                and (inner_product <= 5)
        )

        conditiona3 = (
                lambda_0_e7[0] == lambda_0_e7[1] == lambda_0_e7[2] == lambda_0_e7[3] < lambda_0_e7[4]
                and (inner_product <= 5) and np.allclose(lambda_0_e7[0] % 1, 0)
        )

        conditiona4 = (
                0 < lambda_0_e7[0] == lambda_0_e7[1] == lambda_0_e7[2] < lambda_0_e7[3]
                and (inner_product <= 4) and np.allclose(lambda_0_e7[0] % 1, 0.5)
        )

        conditiona5 = (
                lambda_0_e7[0] == lambda_0_e7[1] < lambda_0_e7[2]
                and (inner_product <= 3)
        )

        conditiona6 = (
                lambda_0_e7[0] - lambda_0_e7[1] - lambda_0_e7[2] - lambda_0_e7[3] - lambda_0_e7[4] - lambda_0_e7[5] -
                lambda_0_e7[6] + lambda_0_e7[7] == 0
                and (inner_product <= 2)
        )

        conditiona7 = (
                lambda_0_e7[0] - lambda_0_e7[1] - lambda_0_e7[2] - lambda_0_e7[3] - lambda_0_e7[4] - lambda_0_e7[
            5] - lambda_0_e7[6] + lambda_0_e7[7] > 0
                and (inner_product <= 1)
        )

        condition2 = (
                lambda_0_e7[0] == lambda_0_e7[1] == lambda_0_e7[2] == lambda_0_e7[3] == lambda_0_e7[4] == 0 and
                lambda_0_e7[5] == -17
                and (inner_product <= 9 or inner_product == 13 or inner_product == 17)
        )

        condition3 = (
                lambda_0_e7[0] == lambda_0_e7[1] == lambda_0_e7[2] == lambda_0_e7[3] == 0
                and (inner_product <= 5 or inner_product == 9) and (lambda_0_e7[4] + lambda_0_e7[5] == -17)
        )

        if (conditiona1 or conditiona2 or conditiona3 or conditiona4 or conditiona5 or conditiona6 or conditiona7 or
                condition2 or condition3):
            filtered_results.append((1, result, inner_product))
            #print(f"-wrho is { result}")
            print(f"-wrho-rho = {(result-rho_e7).tolist()}")
            #print(f"lambda_0 is {lambda_0_e7}")
            print(f"L_w is unitary")
        else:
           # print(f"-wrho is {result}")
           print(f"-wrho-rho = {(result-rho_e7).tolist()}")
           #print(f"lambda_0 is {lambda_0_e7}")
           print(f"L_w is not unitary")
    elif k >= 4:
        print(f"Error")
    else:
        lambda_0_e7 = result - rho_e7 - inner_product * xi_e7
        print(f"-wrho-rho = {(result-rho_e7).tolist()}")
        #print(f"lambda_0 is {lambda_0_e7}")
        print(f"L_w is not unitary")



