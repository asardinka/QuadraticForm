import numpy as np

# Если матрица квадратная
size = 5

rows = n
cols = n
file_name = f"test_{n}.txt"
symmetric = True 

matrix = np.random.uniform(-50000, 50000, size=(rows, cols))
# Округление до двух знаков после запятой
matrix_rounded = np.round(matrix, 2)

if symmetric:
    for i in range(rows):
        for j in range(i + 1, cols):
            matrix_rounded[j, i] = matrix_rounded[i, j]

with open(file_name, 'w') as f:
    for row in matrix_rounded:
        f.write(' '.join(map(str, row)) + '\n')

print(f"Матрица сохранена в '{file_name}'")