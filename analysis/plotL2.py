import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

x = sp.symbols('x')

def lagrange_interpolating_polynomial(x_values, y_values):
    n = len(x_values)
    polynomial = 0
    for i in range(n):

        basis_polynomial = 1
        for j in range(n):
            if i != j:
                basis_polynomial *= (x - x_values[j]) / (x_values[i] - x_values[j])

        polynomial += y_values[i] * basis_polynomial
    return sp.simplify(polynomial)


dx1, dy1 = [9.4000000000000006e-01, 9.4310102051443367e-01, 9.5289897948556646e-01, 9.6000000000000008e-01], [4.3782022467408606e+00, 4.5480278278701407e+00, 5.0000000498219004e+00,  5.0000000497417698e+00]
dx2 = [9.4000000000000006e-01, 9.4155051025721692e-01, 9.4644948974278320e-01, 9.5000000000000007e-01, 9.5155051025721682e-01, 9.5644948974278321e-01, 9.5999999999999996e-01]
dy2 = [4.3802766874711221e+00,4.4626544343208137e+00,4.7542265496785259e+00,4.9999997775742280e+00,5.0000000491017351e+00,5.0000000498325416e+00,5.0000000494853518e+00]

p1 = lagrange_interpolating_polynomial(dx1, dy1)
p2 = lagrange_interpolating_polynomial(dx2[:4], dy2[:4])
p3 = lagrange_interpolating_polynomial(dx2[3:], dy2[3:])

def f1(xval):
    return [p1.subs(x, v) for v in xval]

def f2(xval):
    return [p2.subs(x, v) for v in xval]

def f3(xval):
    return [p3.subs(x, v) for v in xval]


x1 = np.linspace(0.94, 0.96, 400)
x2 = np.linspace(0.94, 0.95, 400)
x3 = np.linspace(0.95, 0.96, 400)

y1 = f1(x1)
y2 = f2(x2)
y3 = f3(x3)


plt.figure(figsize=(10, 6))
plt.scatter(dx1, dy1, color='red', label='Iteration 0')
plt.plot(x1, y1, color='orange', label='Interval (0.92 to 0.94)')

plt.scatter(dx2, dy2, color='blue', label='Iteration 1')
plt.plot(x2, y2, color='green', label='Subinterval 1 (0.94 to 0.95)')
plt.plot(x3, y3, color='lime', label='Subinterval 2 (0.95 to 0.96)')

plt.xlabel('x')
plt.ylabel('y')
plt.title('L2 Coefficient = 0.58')
plt.legend()
plt.grid(True)

plt.show()