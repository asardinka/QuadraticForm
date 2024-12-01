#include "Rational.h"
#include "Matrix.h"
#include "Vector.h"
#include "QuadraticForm.h"
#include "Tests.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>


using namespace std;

int main() {

    // test_eigenvalues();
    // test_eigenvectors();
    //test_Lagrange();
    // test_Lagransh();
    // test_orthogonal_transformations();
    //test_invariants();

/*
    int n = 3;
    QuadraticForm M1(n);
    cin >> M1;
    M1.toSymmetric();
    cout << M1 << endl;

    Vector eigenvalues = M1.getEigenvalues();
    QuadraticForm eigenvectors = M1.getEigenvectors();
    cout <<"eigenvalues: " << endl << eigenvalues << endl;
    cout <<"eigenvectors: " << endl << eigenvectors << endl;

    QuadraticForm M2(n);
    QuadraticForm M3(n);


    //auto start = std::chrono::high_resolution_clock::now();

    M2 = M1.toNormalForm();
    M3 = M1.lagransh();
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed = end - start;
    // std::cout << "Size: "<< n << "x" << n << " | Time: " << elapsed.count()  << std::endl;

    cout << "Normal form: "<< endl << M2 << endl << endl;
    cout << "lagransh: "<< endl << M3 << endl << endl;

*/

    //
    // QuadraticForm test(n);
    // ifstream file("test_matrix.txt");
    // file >> test;
    // file.close();
    //
    // test.toSymmetric();
    //
    //
    // Vector eigenvalues = test.getEigenvalues();
    // QuadraticForm eigenvectors = test.getEigenvectors();
    // cout <<"eigenvalues: " << endl << eigenvalues << endl;
    // cout <<"eigenvectors: " << endl << eigenvectors << endl;
    //
    //
    //
    // QuadraticForm canonicalForm = test.toNormalForm();
    // cout << canonicalForm;


    int n = 10;
    QuadraticForm test(n);
    ifstream file("test_matrix.txt");
    file >> test;
    file.close();
    test.toSymmetric();
    cout << test;
    QuadraticForm canonicalForm = test.toCanonicalForm();
    ofstream f1("test_10.txt");
    ofstream f2("result_10.txt");
    f1 << test;
    f1.close();
    f2 << canonicalForm;






    return 0;
}