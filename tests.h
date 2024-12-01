//
// Created by User on 09.11.2024.
//

#ifndef TESTS_H
#define TESTS_H

#include "QuadraticForm.h"
#include "Vector.h"
#include "Matrix.h"

void test_eigenvalues() {
    for (int n = 100; n <= 500; n+=100) {
        QuadraticForm test(n);
        ifstream file("test_matrix.txt");
        file >> test;
        file.close();
        test.toSymmetric();

        auto start = chrono::high_resolution_clock::now();

        Vector eigenvalues = test.getEigenvalues();

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end - start;
        cout << "Testing eigenvalues: " << endl;
        cout << "Size: "<< n << "x" << n << " | Time: " << elapsed.count()  << endl << endl;
    }
}
void test_eigenvectors() {
    for (int n = 100; n <= 500; n+=100) {
        QuadraticForm test(n);
        ifstream file("test_matrix.txt");
        file >> test;
        file.close();
        test.toSymmetric();

        auto start = chrono::high_resolution_clock::now();

        QuadraticForm eigenvectors = test.getEigenvectors();

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end - start;
        cout << "Testing eigenvectors: " << endl;
        cout << "Size: "<< n << "x" << n << " | Time: " << elapsed.count()  << endl << endl;
    }
}
void test_orthogonal_transformations() {
    for (int n = 100; n <= 500; n+=100) {
        QuadraticForm test(n);
        ifstream file("test_matrix.txt");
        file >> test;
        file.close();
        test.toSymmetric();

        auto start = chrono::high_resolution_clock::now();

        QuadraticForm canonicalForm = test.toCanonicalForm();

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end - start;
        cout << "Testing method of orthogonal transformations: " << endl;
        cout << "Size: "<< n << "x" << n << " | Time: " << elapsed.count()  << endl << endl;
    }
}
void test_Lagrange() {
    for (int n = 500; n <= 5000; n+=500) {
        QuadraticForm test(n);
        ifstream file("test_matrix.txt");
        file >> test;
        file.close();
        test.toSymmetric();

        auto start = chrono::high_resolution_clock::now();

        QuadraticForm canonicalForm = test.toNormalForm();

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end - start;
        cout << "Testing Lagrange method: " << endl;
        cout << "Size: "<< n << "x" << n << " | Time: " << elapsed.count()  << endl << endl;
    }
}
void test_Lagransh() {
    for (int n = 100; n <= 1000; n+=100) {
        QuadraticForm test(n);
        ifstream file("test_matrix.txt");
        file >> test;
        file.close();
        test.toSymmetric();

        auto start = chrono::high_resolution_clock::now();

        QuadraticForm canonicalForm = test.lagransh();

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end - start;
        cout << "Testing Lagransh method: " << endl;
        cout << "Size: "<< n << "x" << n << " | Time: " << elapsed.count()  << endl << endl;
    }
}
void test_invariants() {

    for (int n = 5; n <= 10; n += 5) {
        QuadraticForm test(n);
        ifstream file("test_matrix.txt");
        file >> test;
        file.close();
        test.toSymmetric();
        cout << "Testing for the preservation of invariants on matrix with size " << n << " x " << n << endl;

        Vector eigenvalues = test.getEigenvalues();

        // First invariants

        double Trace = test.Trace();
        double sumEigenvalues = eigenvalues.sum();
        cout << "Trace: " << Trace << endl;
        cout << "Sum eigenvalues: " << sumEigenvalues << endl << endl;

        // Second invariants

        double sumSquares = test.sum_of_squares();
        double sum_squares_Eigenvalues = eigenvalues.sumOfSquares();
        cout << "Sum of squares: " << sumSquares << endl;
        cout << "Sum of squares eigenvalues: : " << sum_squares_Eigenvalues << endl << endl;

        // Third invariants
        double det = test.determinant();
        double prodEigenvalues = eigenvalues.product();
        cout << "Determinant: " << det << endl;
        cout << "Product eigenvalues: " << prodEigenvalues << endl << endl;

        cout << "____________________________________" << endl;
    }
}

#endif //TESTS_H
