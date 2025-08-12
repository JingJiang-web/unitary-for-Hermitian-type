import numpy as np
from fractions import Fraction
Coxeter_type = input("Coxeter type (e.g., A, B, C, D, D^star): ")
n = int(input("n : "))
w_input = input("w (e.g., [1, 3, 4]): ")
try:
    w = [int(i) for i in w_input.strip('[]').split(',') if i.strip()]
except ValueError:
    print("The form of w is incorrect")
    exit(1)

if Coxeter_type == 'A':
    p = int(input("p (the number of entries in the second column of Young tableau): "))
    rho = np.array([n / 2 - i for i in range(n + 1)])
    beta = np.array([1] + [0] * (n - 1) + [-1])
    xi = np.array([(n + 1 - p) / (n + 1)] * (p)+ [-1 * p / (n + 1)] * (n + 1 - p))
elif Coxeter_type == 'B':
    rho = np.array([n - i - 1 / 2 for i in range(n)])
    beta = np.array([1] + [1] + [0] * (n - 2))
    xi = np.array([1] + [0] * (n - 1))
elif Coxeter_type == 'C':
    k = int(input("k : "))
    rho = np.array([n - i for i in range(n)])
    beta = np.array([2] + [0] * (n - 1))
    xi = np.array([1] * (n))
elif Coxeter_type == 'D':
    rho = np.array([n - i -1 for i in range(n)])
    beta = np.array([1] + [1] + [0] * (n - 2))
    xi = np.array([1] + [0] * (n - 1))
elif Coxeter_type == 'D^star':
    rho = np.array([n - i -1 for i in range(n)])
    beta = np.array([1] + [1] + [0] * (n - 2))
    xi = np.array([1 / 2] * (n))
else:
    print("Coxeter type is invalid")
    exit(1)

# 构建矩阵
def build_matrix(n, last_diag = 1, extra_col = False):
    if extra_col:
        matrix = [[0] * (n + 1) for _ in range(n)]
    else:
        matrix = [[0] * n for _ in range(n)]
    for i in range(n):
        matrix[i][i] = 1
        if i < n - 1:
            matrix[i][i + 1] = -1
    if extra_col:
        if Coxeter_type == 'A':
            matrix[n - 1][n] = -1
        elif Coxeter_type == 'D' or 'D^star':
            matrix[n - 1] = [0] * (n - 2) + [1, 1] + ([1] if extra_col else [])
            matrix[n - 1][n] = 1
    if last_diag != 1:
        matrix[n - 1][n - 1] = last_diag
    return matrix

# 反射函数
def reflection(x, a):
    if len(x) != len(a):
        print("Error")
        return x
    dot_a_a = np.dot(a, a)
    if dot_a_a == 0:
        return x
    return x - 2 * np.dot(a, x) / dot_a_a * a

if Coxeter_type == 'A':
    matrix = build_matrix(n, extra_col = True)
    a = [np.array(matrix[i])[:n+1] for i in range(n)]
elif Coxeter_type == 'B':
    matrix = build_matrix(n)
    a = [np.array(matrix[i]) for i in range(n)]
elif Coxeter_type == 'C':
    matrix = build_matrix(n, last_diag = 2)
    a = [np.array(matrix[i]) for i in range(n)]
elif Coxeter_type == 'D':
    matrix = build_matrix(n, extra_col = True)
    a = [np.array(matrix[i])[:n] for i in range(n)]
elif Coxeter_type == 'D^star':
    matrix = build_matrix(n, extra_col = True)
    a = [np.array(matrix[i])[:n] for i in range(n)]

reflection_functions = [lambda x, a=a: reflection(x, a) for a in a]

def compose(*functions):
    def composer(value):
        for function in reversed(functions):
            value = function(value)
        return value
    return composer

valid_w = [i - 1 for i in w if 1 <= i <= len(reflection_functions)]
if not valid_w:
    print("w is invalid")
    exit(1)

def find_p1_q1(arr):
    if len(arr) == 0:
        return 0, 0
    # 计算 p1
    first_element = arr[0]
    p1 = 0
    for element in arr:
        if element == first_element:
            p1 += 1
        else:
            break
    last_element = arr[-1]
    q1 = 0
    for element in reversed(arr):
        if element == last_element:
            q1 += 1
        else:
            break
    return p1, q1

def find_q_r(arr):
    if len(arr) == 0:
        return 0, 0
    first_element = arr[0]
    q = 0
    for element in arr:
        if element == first_element:
            q = q + 1
        else:
            break
    target_num = first_element - 1
    r = 0
    for idx in range(q, len(arr)):
        if arr[idx] == target_num:
            r = idx + 1
    if r == 0:
        r = q
    return q, r

composed_function = compose(*[reflection_functions[i] for i in valid_w])
result = -1 * composed_function(rho)
inner_product = 2 * np.dot(result, beta) / np.dot(beta, beta)
lambda_0 = result - rho - inner_product * xi

if Coxeter_type == 'A':
    if (inner_product <= n):
        p1, q1 = find_p1_q1(lambda_0)
        #print(p1, q1)
        if (inner_product <= max(p1, q1)) or (isinstance (inner_product, int) and inner_product <= p1 + q1 - 1):
            lambda_0_fractions = [Fraction(str(num)).limit_denominator() for num in lambda_0]
            str_arr = []
            for fraction in lambda_0_fractions:
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
            print(f"-wrho-rho = {(result - rho).tolist()}")
            print(f"-wrho = {(result).tolist()}")
            print(f"inner product is {inner_product}")
            print(f"lambda_0 is {str_arr}")
            print(f"L_w is unitary")
            print("-----------------------------------------------------------------")
        else:
            lambda_0_fractions = [Fraction(str(num)).limit_denominator() for num in lambda_0]
            str_arr = []
            for fraction in lambda_0_fractions:
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
            print(f"-wrho-rho = {(result - rho).tolist()}")
            print(f"-wrho = {(result).tolist()}")
            print(f"inner product is {inner_product}")
            print(f"lambda_0 is {str_arr}")
            print(f"L_w is not unitary")
            print("-----------------------------------------------------------------")

elif Coxeter_type == 'B':
    condition2 = (
        all(lambda_0[j] == 0 for j in range(1, n))
        and (inner_product <= n - 1 / 2 or inner_product == 2 * n - 2)
    )
    condition3 = (
        all(lambda_0[j] == 1 / 2 for j in range(1, n))
        and inner_product <= n - 1 / 2
    )
    if any(lambda_0[1] == lambda_0[p] > lambda_0[p+1] > 0 and inner_product <= p for p in range(1, n)) or condition2 or condition3:
        print(f"-wrho-rho = {(result - rho).tolist()}")
        print(f"-wrho = {(result).tolist()}")
        print(f"inner product is {inner_product}")
        print(f"lambda_0 = {(lambda_0).tolist()}")
        print(f"L_w is unitary")
        print("-----------------------------------------------------------------")
    else:
        print(f"-wrho-rho = {(result - rho).tolist()}")
        print(f"-wrho = {(result).tolist()}")
        print(f"inner product is {inner_product}")
        print(f"lambda_0 = {(lambda_0).tolist()}")
        print(f"L_w is not unitary")
        print("-----------------------------------------------------------------")

elif Coxeter_type == 'C':
    q, r = find_q_r(lambda_0)
    if (k < n and inner_product == n - k / 2 ) and (inner_product <= (r + 1) / 2 or (isinstance (2 * inner_product, int) and inner_product <= (q + r) / 2)):
        print(f"-wrho-rho = {(result - rho).tolist()}")
        print(f"-wrho = {(result).tolist()}")
        print(f"inner product is {inner_product}")
        print(f"lambda_0 = {(lambda_0).tolist()}")
        print(f"L_w is unitary")
        print("-----------------------------------------------------------------")
    elif (k == n) and (inner_product <= (r + 1) / 2 or (isinstance (2 * inner_product, int) and inner_product <= (q + r) / 2)):
        print(f"-wrho-rho = {(result - rho).tolist()}")
        print(f"-wrho = {(result).tolist()}")
        print(f"inner product is {inner_product}")
        print(f"lambda_0 = {(lambda_0).tolist()}")
        print(f"L_w is unitary")
        print("-----------------------------------------------------------------")

    else:
        print(f"-wrho-rho = {(result - rho).tolist()}")
        print(f"-wrho = {(result).tolist()}")
        print(f"inner product is {inner_product}")
        print(f"lambda_0 = {(lambda_0).tolist()}")
        print(f"L_w is not unitary")
        print("-----------------------------------------------------------------")

elif Coxeter_type == 'D^star':
    q, r = find_q_r(lambda_0)
    if (any(lambda_0[0] >= lambda_0[1] == lambda_0[p] >= lambda_0[p + 1] and inner_product <= p for p in range(1, n))) or (q % 2 == 1 and q >= 3 and inner_product <= 2 * q - 3 and inner_product != q + 1) or (q % 2 == 0 and q >= 3 and inner_product <= 2 * q - 3 and inner_product != q):
        #print(q)
        print(f"-wrho-rho = {(result - rho).tolist()}")
        print(f"-wrho = {(result).tolist()}")
        print(f"inner product is {inner_product}")
        print(f"lambda_0 = {(lambda_0).tolist()}")
        print(f"L_w is unitary")
        print("-----------------------------------------------------------------")
    else:
        print(f"-wrho-rho = {(result - rho).tolist()}")
        print(f"-wrho = {(result).tolist()}")
        print(f"inner product is {inner_product}")
        print(f"lambda_0 = {(lambda_0).tolist()}")
        print(f"L_w is not unitary")
        print("-----------------------------------------------------------------")

elif Coxeter_type == 'D':
    if (any(lambda_0[1] == lambda_0[p] > lambda_0[p+1] >= 0 and inner_product <= p for p in range(1, n))) or (lambda_0[1] == lambda_0[n-1] == 0 and (inner_product <= n - 1 or inner_product == 2 * n - 3)):
        print(f"-wrho-rho = {(result - rho).tolist()}")
        print(f"-wrho = {(result).tolist()}")
        print(f"inner product is {inner_product}")
        print(f"lambda_0 = {(lambda_0).tolist()}")
        print(f"L_w is unitary")
        print("-----------------------------------------------------------------")
    else:
        print(f"-wrho-rho = {(result - rho).tolist()}")
        print(f"-wrho = {(result).tolist()}")
        print(f"inner product is {inner_product}")
        print(f"lambda_0 = {(lambda_0).tolist()}")
        print(f"L_w is not unitary")
        print("-----------------------------------------------------------------")