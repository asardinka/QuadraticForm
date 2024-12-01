//
// Created by Admin on 19.10.2024.
//

#include "QuadraticForm.h"
#include  <cmath>

QuadraticForm::QuadraticForm(const Matrix& M){
    row = M.getRow();
    col = M.getCol();
    isTranspose = false;
    ptr = new double*[row];
    for (int i = 0; i < row; i++) {
        ptr[i] = new double[col];
        for (int j = 0; j < col; j++)
            ptr[i][j] = M(i,j);
    }
}

QuadraticForm::QuadraticForm(int size) {
    row = size;
    col = size;
    isTranspose = false;
    ptr = new double*[row];
    for (int i = 0; i < row; ++i)
        ptr[i] = new double[col];
}

int QuadraticForm::getSize() const {return row;}

QuadraticForm QuadraticForm::gramSchmidtOrthogonalization() {
    //ptr[0][0] += 0.0000001;
    QuadraticForm orthogonalMatrix(row);
    for (int i = 0; i < col; i++) {
        Vector currentVector(row);
        for (int j = 0; j < row; j++)
            currentVector[j] = (*this)(j, i);
        for (int k = 0; k < i; k++) {
            Vector previousVector(row);
            for (int j = 0; j < row; j++)
                previousVector[j] = orthogonalMatrix(j, k);
            double projection = currentVector.dot(previousVector) / previousVector.dot(previousVector);
            for (int j = 0; j < row; j++)
                currentVector[j] = currentVector[j] - (projection * previousVector[j]);
        }
        double norm = currentVector.norm();
        if (norm > double(0)) {
            for (int j = 0; j < row; j++)
                orthogonalMatrix(j, i) = currentVector[j] / norm;
        } else {
            for (int j = 0; j < row; j++)
                orthogonalMatrix(j, i) = 0;
        }
    }
    return orthogonalMatrix;
}


void QuadraticForm::QRdecomposition(Matrix& Q, Matrix& R) {
    QuadraticForm q = gramSchmidtOrthogonalization();
    Q = q;
    q.transpose();
    R = q * (*this);
}

Vector QuadraticForm::getEigenvalues() {
    Vector eigenvalues(row);
    QuadraticForm A = *this;
    int maxIterations = 100;
    int iteration = 0;
    for (; iteration < maxIterations; ++iteration) {
        Matrix Q(row, col);
        Matrix R(row, col);
        A.QRdecomposition(Q, R);
        A = R * Q;
        if (A.isTriangular()) {
            break;
        }
    }
    for (int i = 0; i < row; ++i)
        eigenvalues[i] = A(i, i);
    return eigenvalues;
}

QuadraticForm QuadraticForm::getEigenvectors() {
    QuadraticForm result(row);
    QuadraticForm A = *this;
    Matrix Q(row, row);
    Matrix R(row, row);
    for (int i = 0; i < row; ++i){
        for (int j = 0; j < col; ++j) {
            if(i == j)
                result(i,j) = 1;
            else
                result(i, j) = 0;
        }
    }
    int maxIterations = 50;
    for (int iteration = 0; iteration < maxIterations; ++iteration) {
        A.QRdecomposition(Q, R);
        A = R * Q;
        result = result * Q;
        if (A.isTriangular())
            break;
    }

    for (int i = 0; i < col; i++){
        double c = result(0,i);
        for (int j = 0; j < row; j++)
            result(j,i) /= abs(c);
    }

    return result;
}

QuadraticForm QuadraticForm::toCanonicalForm(){
    QuadraticForm canonicalForm(row);
    Vector eigenvalues = getEigenvalues();
    for (int i = 0; i < row; i++) {
        double roundedEigenvalue =  round(eigenvalues[i] * double(1000)) / double(1000);
        canonicalForm(i, i) = roundedEigenvalue;
        for (int j = 0; j < col; j++) {
            if (i != j)
                canonicalForm(i, j) = double(0);
        }
    }
    return canonicalForm;
}

QuadraticForm QuadraticForm::toNormalForm() {
    Matrix canonicalForm = *this;
    for (int i = 0; i < row - 1; i++) {
        if (canonicalForm(i, i) == 0) {
            bool foundNonZero = false;
            for (int k = i + 1; k < row; k++) {
                if (canonicalForm(k, i) != 0) {
                    canonicalForm.swapRows(i, k);
                    canonicalForm.swapCols(i, k);
                    foundNonZero = true;
                    break;
                }
            }
            if (!foundNonZero) {
                continue;
            }
        }
        bool hasSquareTerms = false;
        for (int j = i; j < row; j++) {
            if (canonicalForm(j, j) != 0) {
                hasSquareTerms = true;
                break;
            }
        }
        if (!hasSquareTerms) {
            bool isTransform = false;
            for (int k = i; k < row - 1 && !isTransform; k++) {
                for (int l = k + 1; l < row; l++) {
                    if (canonicalForm(k, l) != 0) {
                        // Специальное невырожденное линейное преобразование
                        double a = sqrt(abs(canonicalForm(k, l)));
                        for (int m = 0; m < row; m++) {
                            canonicalForm(k, m) += a * canonicalForm(l, m);
                            canonicalForm(l, m) -= a * canonicalForm(k, m);
                        }
                        isTransform = true;
                        break;
                    }
                }
            }
        }
        for (int k = i + 1; k < row; k++) {
            if (canonicalForm(i, i) == 0) continue;

            double c = canonicalForm(k, i) / canonicalForm(i, i);
            for (int j = i; j < col; j++) {
                canonicalForm(k, j) -= c * canonicalForm(i, j);
            }
        }
        for (int j = 0; j < col; j++) {
            if (i == j && canonicalForm(i, j) != 0)
                canonicalForm(i, j) /= abs(canonicalForm(i, j));
            else
                canonicalForm(i, j) = 0;
        }
    }
    if (canonicalForm(row - 1, col - 1) != 0)
        canonicalForm(row - 1, col - 1) /= abs(canonicalForm(row - 1, col - 1));
    for (int j = 0; j < col - 1; j++) {
        canonicalForm(row - 1, j) = 0;
    }

    return canonicalForm;
}

QuadraticForm QuadraticForm::lagransh(ostream* stream) {
    QuadraticForm currentMatrix = *this;
    QuadraticForm copyCurrent = currentMatrix;
    int n = this->row;
    QuadraticForm result(n);
    for(int x=0;x<n;x++)
        for(int y=0;y<n;y++)
            result(x,y) = 0;
    for(int i=0;i<n;i++){
        copyCurrent = currentMatrix;
        result(i,i) = copyCurrent(i,i);
        for(int x=0;x<n-i;x++)
            for(int y=0;y<n-i;y++)
                currentMatrix(x+i,y+i) -= copyCurrent(x+i,i)*copyCurrent(i,y+i)*(1/copyCurrent(i,i));
    }
    //    for(int x=0;x<n;x++)
    //        for(int y=0;y<n;y++)
    //            if (x!=y)
    //                result(x,y) = 0;
    return result;
}
double QuadraticForm::sum_of_squares() const {
    double sum = 0.0;
    for (int i = 0; i < row; ++i) {
        sum += (*this)(i, i) * (*this)(i, i);
    }
    for (int i = 0; i < row; ++i) {
        for (int j = i + 1; j < col; ++j) {
            sum += 2 * (*this)(i, j) * (*this)(i, j);
        }
    }
    return sum;
}



