import numpy as np

# Define test matrices
test_matrices = {
    "Identity Matrix": np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
    "Diagonal Matrix": np.array([[3, 0, 0], [0, 5, 0], [0, 0, 7]]),
    "Symmetric Matrix": np.array([[6, 2, 0], [2, 3, 0], [0, 0, 5]]),
    "Random Symmetric Matrix": np.array([[4, -1, 0.5], [-1, 2, 0], [0.5, 0, 3]]),
}

# Perform eigen decomposition using NumPy
for name, matrix in test_matrices.items():
    print(f"Test: {name}")
    print("Matrix:")
    print(matrix)

    eigenvalues, eigenvectors = np.linalg.eig(matrix)
    print("Eigenvalues:")
    print(eigenvalues)
    print("Eigenvectors:")
    print(eigenvectors)
    print("-" * 50)
