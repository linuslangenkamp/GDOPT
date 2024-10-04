### basically chatgpt code

from sympy import symbols, Matrix, eye
import numpy as np

import matplotlib.pyplot as plt

"""
c | A
-----
    b
"""

# Define z as a symbolic variable
z_sym = symbols("z")


# custom 2 step order 2
# Redefine the matrix A using symbolic z
"""
4/7 | G   :  G-4/7
 1  | 7/6 :  -1/6
------------------
	| 7/6 :  -1/6
"""
"""
G = 4/7
A_sym = Matrix([
    [G, 4/7-G],
    [7/6, -1/6],
])	

# Define the b vector and ones vector symbolically
b_sym = Matrix([7/6, -1/6])


-> general c = [n / (2n + 1), 1] for some n
-> a21, a22 != b -> b = a2: = [(2n+1)/(2n+2), 1/(2n+2)];; RadauIIA 2 step 1/3, 1 -> 3/4 1/4
"""
"""
RK 4

A_sym = Matrix([
    [0, 0, 0, 0],
    [0.5, 0, 0, 0],
    [0, 0.5, 0, 0],
    [0, 0, 1, 0]
])

# Define the b vector using floats
b_sym = Matrix([
    1/6,
    1/3,
    1/3,
    1/6
])
"""


c1 = 0.5
det = 0.5 / (c1 - 1)
A_sym = Matrix(
    [
        [(-2 * c1 + c1 * c1) / det, c1 * c1 / det],
        [-1.0 / det, (2.0 * c1 - 1.0) / det],
    ]
)

# Define the b vector and ones vector symbolically
b_sym = Matrix([-1.0 / det, (2.0 * c1 - 1.0) / det])


# Define the identity matrix I
I_sym = eye(len(b_sym))

# Compute I - zA symbolically
IzA_sym = I_sym - z_sym * A_sym

# Invert the matrix I - zA
IzA_inv_sym = IzA_sym.inv()

ones_sym = Matrix([1] * len(b_sym))

# Compute the stability function symbolically: R(z) = 1 + z * b.T * (I - zA)^(-1) * ones
R_z_sym = 1 + z_sym * (b_sym.T * (IzA_inv_sym * ones_sym))[0]

R_z_sym.simplify()

print(R_z_sym.simplify())

from sympy import lambdify

# Convert the symbolic R(z) to a numerical function that works with NumPy arrays
R = lambdify(z_sym, R_z_sym, "numpy")

# Create a grid of complex numbers
real_vals = np.linspace(-10, 10, 250)
imag_vals = np.linspace(-10, 10, 250)
Z_real, Z_imag = np.meshgrid(real_vals, imag_vals)
Z = Z_real + 1j * Z_imag

# Compute |R(z)| for each point in the grid
R_values = R(Z)
abs_R = np.abs(R_values)

# Create a contour plot for |R(z)| <= 1
plt.figure(figsize=(10, 8))

# Plot only the region where |R(z)| <= 1 using the masked array
plt.contourf(Z_real, Z_imag, abs_R, levels=np.linspace(0, 1, 500), cmap="RdYlBu", vmin=0, vmax=1)

# Add a colorbar for the plot with correct labels
plt.colorbar(label=r"$|R(z)| \, \text{for} \, |R(z)| \leq 1$")

# Add a contour line for |R(z)| = 1 (the boundary of stability)
plt.contour(Z_real, Z_imag, abs_R, levels=[1], colors="black", linewidths=1.5)

# Plot settings and labels
plt.title(r"Region of Stability ($|R(z)| \leq 1$)")
plt.xlabel("Real part of $z$")
plt.ylabel("Imaginary part of $z$")
plt.grid()
plt.xlim(-10, 10)
plt.ylim(-10, 10)
plt.axhline(0, color="black", linewidth=0.5, linestyle="--")
plt.axvline(0, color="black", linewidth=0.5, linestyle="--")

plt.show()
