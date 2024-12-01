import numpy as np

rows = 5
cols = 5
file_name = "test_matrix.txt"
symmetric = True 

matrix = np.random.randint(-1000, 1000, size=(rows, cols))

if symmetric:
    for i in range(rows):
        for j in range(i + 1, cols):
            matrix[j, i] = matrix[i, j]

with open(file_name, 'w') as f:
    for row in matrix:
        f.write(' '.join(map(str, row)) + '\n')

print(f"Матрица сохранена в '{file_name}'")