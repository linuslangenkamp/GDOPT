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
z = symbols("z")

A = Matrix(
    [
        [5/12, -1/12],
        [3/4, 1/4],
    ]
)

# Define the b vector and ones vector symbolically
b = Matrix([3/4, 1/4])

size = 5


inverse = (eye(len(b)) - z * A).inv()

ones = Matrix([1] * len(b))

# Compute the stability function symbolically: R(z) = 1 + z * b.T * (I - zA)^(-1) * ones
Rz = 1 + z * (b.T * (inverse * ones))[0]

Rz.simplify()

print(Rz.simplify())

from sympy import lambdify

# Convert the symbolic R(z) to a numerical function that works with NumPy arrays
R = lambdify(z, Rz, "numpy")

# Create a grid of complex numbers
realVals = np.linspace(-size, size, 500)
imagVals = np.linspace(-size, size, 500)
Zreal, Zimag = np.meshgrid(realVals, imagVals)
Z = Zreal + 1j * Zimag

# Compute |R(z)| for each point in the grid
Rvals = R(Z)
Rabs = np.abs(Rvals)

# Create a contour plot for |R(z)| <= 1
plt.figure(figsize=(10, 8))

# Plot only the region where |R(z)| <= 1 using the masked array
# add cmap
plt.contourf(Zreal, Zimag, Rabs, levels=np.linspace(0, 1, 500), cmap="coolwarm", vmin=0, vmax=1)

# Add a colorbar for the plot with correct labels
plt.colorbar(label=r"$|R(z)| \, \text{for} \, |R(z)| \leq 1$")

# Add a contour line for |R(z)| = 1 (the boundary of stability)
plt.contour(Zreal, Zimag, Rabs, levels=[1], colors="black", linewidths=1.5)

# Plot settings and labels
plt.title(r"Region of Stability ($|R(z)| \leq 1$)")
plt.xlabel("Real part of $z$")
plt.ylabel("Imaginary part of $z$")
plt.grid()
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.axhline(0, color="black", linewidth=0.5, linestyle="--")
plt.axvline(0, color="black", linewidth=0.5, linestyle="--")

plt.show()
