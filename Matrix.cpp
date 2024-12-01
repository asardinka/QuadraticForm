#include "Matrix.h"
#include  "math.h"
Matrix::Matrix(int r, int c) : row(r), col(c), isTranspose(false) {
    ptr = new double*[row];
    for (int i = 0; i < row; ++i) {
        ptr[i] = new double[col];
    }
}

Matrix::Matrix(const Matrix& M) : row(M.row), col(M.col), isTranspose(M.isTranspose) {
    ptr = new double*[row];
    for (int i = 0; i < row; i++) {
        ptr[i] = new double[col];
        for (int j = 0; j < col; j++) {
            ptr[i][j] = M.ptr[i][j];
        }
    }
}

Matrix::~Matrix() {
    if (ptr != nullptr) {
        for (int i = 0; i < row; i++) {
            delete[] ptr[i];
        }
        delete[] ptr;
        ptr = nullptr;
    }
}

void Matrix::transpose() {
    isTranspose = !isTranspose;
}

double& Matrix::operator()(int i, int j) {
    return isTranspose ? ptr[j][i] : ptr[i][j];
}

const double& Matrix::operator()(int i, int j) const {
    return isTranspose ? ptr[j][i] : ptr[i][j];
}

Matrix& Matrix::operator=(const Matrix& M) {
    if (this == &M) return *this;
    if (row != M.row || col != M.col) {
        for (int i = 0; i < row; i++) {
            delete[] ptr[i];
        }
        delete[] ptr;
        row = M.row;
        col = M.col;
        ptr = new double*[row];
        for (int i = 0; i < row; i++) {
            ptr[i] = new double[col];
        }
    }
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            ptr[i][j] = M.ptr[i][j];
        }
    }
    return *this;
}

double Matrix::Trace() {
    double result = 0;
    for (int i = 0; i < std::min(row, col); ++i) {
        result = result + (*this)(i, i);
    }
    return result;
}

Matrix Matrix::operator+(const Matrix& M) {
    Matrix result(row, col);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            result(i, j) = (*this)(i, j) + M(i, j);
        }
    }
    return result;
}

Matrix Matrix::operator-(const Matrix& M) {
    Matrix result(row, col);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            result(i, j) = (*this)(i, j) - M(i, j);
        }
    }
    return result;
}

Matrix Matrix::operator*(const double& n) {
    Matrix result(row, col);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            result(i, j) = (*this)(i, j) * n;
        }
    }
    return result;
}

Matrix Matrix::operator*(const Matrix& M) {
    if (col != M.row) {
        throw std::invalid_argument("The sizes of the matrices do not match for multiplication.");
    }
    Matrix result(row, M.col);
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < M.col; ++j) {
            result(i, j) = 0 ;
            for (int k = 0; k < col; ++k) {
                result(i, j) = result(i, j) + ((*this)(i, k) * M(k, j));
            }
        }
    }
    return result;
}

bool Matrix::isTriangular(double epsilon) const {
    int n = row;
    for (int i = 1; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            if (abs((*this)(i, j)) > epsilon) {
                return false;
            }
        }
    }
    return true;
}

bool Matrix::isSymmetric() const {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            if (((*this)(i, j) - (*this)(j, i)) < 0.000001) {
                return false;
            }
        }
    }
    return true;
}

void Matrix::toSymmetric() {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            (*this)(i, j) = (*this)(j, i) = ((*this)(i, j) + (*this)(j, i)) / 2;
        }
    }
}

void Matrix::swapRows(int row1, int row2) {
    if (row1 < 0 || row1 >= row || row2 < 0 || row2 >= row) {
        throw std::out_of_range("Invalid row index for swapRows");
    }
    for (int j = 0; j < col; j++) {
        std::swap((*this)(row1, j), (*this)(row2, j));
    }
}

void Matrix::swapCols(int col1, int col2) {
    if (col1 < 0 || col1 >= col || col2 < 0 || col2 >= col) {
        throw std::out_of_range("Invalid column index for swapCols");
    }
    for (int i = 0; i < row; i++) {
        std::swap((*this)(i, col1), (*this)(i, col2));
    }
}

std::ostream& operator<<(std::ostream& s, Matrix& M) {
    for (int i = 0; i < M.row; i++) {
        for (int j = 0; j < M.col; j++) {
            s << M(i, j) << " ";
        }
        s << std::endl;
    }
    return s;
}
double Matrix::determinant() {
    if (row != col) {
        throw std::invalid_argument("Matrix must be square to calculate determinant.");
    }
    Matrix A(*this);
    double det = 1.0;
    for (int i = 0; i < A.row; ++i) {
        int maxRow = i;
        for (int j = i + 1; j < A.row; ++j) {
            if (abs(A(j, i)) > abs(A(maxRow, i))) {
                maxRow = j;
            }
        }
        if (abs(A(maxRow, i)) < 1e-9) {
            return 0.0;
        }
        if (maxRow != i) {
            A.swapRows(i, maxRow);
            det = -det;
        }
        for (int j = i + 1; j < A.row; ++j) {
            if (abs(A(j, i)) > 1e-9) {
                double factor = A(j, i) / A(i, i);
                for (int k = i; k < A.col; ++k) {
                    A(j, k) -= factor * A(i, k);
                }
            }
        }
    }
    for (int i = 0; i < A.row; ++i) {
        det *= A(i, i);
    }

    return det;
}
std::istream& operator>>(std::istream& s, Matrix& M) {
    for (int i = 0; i < M.row; i++) {
        for (int j = 0; j < M.col; j++) {
            s >> M(i, j);
        }
    }
    return s;
}

std::ifstream& operator>>(std::ifstream& s, Matrix& M) {
    for (int i = 0; i < M.row; i++) {
        for (int j = 0; j < M.col; j++) {
            s >> M(i, j);
        }
    }
    return s;
}

bool Matrix::getTranspose() const {
    return isTranspose;
}

int Matrix::getRow() const {
    return row;
}

int Matrix::getCol() const {
    return col;
}