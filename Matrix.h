//
// Created by Admin on 18.10.2024.
//

#ifndef KURSOVAYA_MATRIX_H
#define KURSOVAYA_MATRIX_H

#include <iostream>
#include <chrono>
#include <fstream>
#include "Vector.h"

class Matrix {
protected:
    double** ptr;
    int col;
    bool isTranspose;
    int row;
public:
    Matrix(int r = 0, int c = 0);
    Matrix(const Matrix& M);
    ~Matrix();

    int getRow() const;
    int getCol() const;
    bool getTranspose() const;

    double& operator()(int i, int j);
    const double& operator()(int i, int j) const;
    Matrix& operator=(const Matrix& M);
    Matrix operator+(const Matrix& M);
    Matrix operator-(const Matrix& M);
    Matrix operator*(const double& n);
    Matrix operator*(const Matrix& M);

    double determinant();
    double Trace();
    void transpose();
    bool isSymmetric() const;
    void toSymmetric();
    bool isTriangular(double epsilon = 1e-8) const;
    void swapRows(int row1, int row2);
    void swapCols(int col1, int col2);

    friend std::ostream& operator<<(std::ostream& s, Matrix& M);
    friend std::istream& operator>>(std::istream& s, Matrix& M);
    friend std::ifstream& operator>>(std::ifstream& s, Matrix& M);
};

#endif //KURSOVAYA_MATRIX_H
