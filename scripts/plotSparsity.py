import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import os

def plotSparsity(m, ax, file_path):
    if not isinstance(m, sp.coo_matrix):
        m = sp.coo_matrix(m)
    for (x, y, data) in zip(m.col, m.row, m.data):
        ax.add_patch(Rectangle(
            xy=(x, y), width=1, height=1, edgecolor='black', facecolor='blue', alpha=0.6))
        """
        if 'hessian' in os.path.basename(file_path).lower() and x != y:
            ax.add_patch(Rectangle(
            xy=(y, x), width=1, height=1, edgecolor='black', facecolor='blue', alpha=0.6)) """
    ax.set_xlim(0, m.shape[1])
    ax.set_ylim(0, m.shape[0])
    ax.invert_yaxis()
    
    if 'jacobian' in os.path.basename(file_path).lower():
        ax.set_xlabel('Variables')
        ax.set_ylabel('Equations')
        ax.set_title('Jacobian Sparsity')
    elif 'hessian' in os.path.basename(file_path).lower():
        ax.set_xlabel('Variables')
        ax.set_ylabel('Variables')
        ax.set_title('Hessian Sparsity')


def plotFromFile(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    dim_row, dim_col = map(int, lines[1].strip().split(','))

    rows, cols = [], []
    for line in lines[3:]:
        row, col = map(int, line.strip().split(','))
        rows.append(row)
        cols.append(col)

    data = np.ones(len(rows))
    sparse_matrix = sp.coo_matrix((data, (rows, cols)), shape=(dim_row, dim_col))

    fig, ax = plt.subplots()
    plotSparsity(sparse_matrix, ax, file_path)
    plt.show()

file_path_j = "/mnt/c/Users/Linus/Desktop/Studium/Master/Masterarbeit/Sparsity/jacobianSparsity.csv"
file_path_h = "/mnt/c/Users/Linus/Desktop/Studium/Master/Masterarbeit/Sparsity/hessianSparsity.csv"
plotFromFile(file_path_h)
plotFromFile(file_path_j)
