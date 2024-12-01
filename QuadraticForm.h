//
// Created by Admin on 19.10.2024.
//

#ifndef KURSOVAYA_QUADRATICFORM_H
#define KURSOVAYA_QUADRATICFORM_H

#include "QuadraticForm.h"
#include "Vector.h"
#include "Matrix.h"

class QuadraticForm : public Matrix {
public:
    QuadraticForm(int size);
    QuadraticForm(const Matrix& matrix);
    int getSize() const;
    QuadraticForm gramSchmidtOrthogonalization();
    void QRdecomposition(Matrix& Q, Matrix& R);
    double sum_of_squares() const;
    Vector getEigenvalues();
    QuadraticForm getEigenvectors();
    QuadraticForm toCanonicalForm();
    QuadraticForm toNormalForm();
    QuadraticForm lagransh(ostream* stream = nullptr);

};


#endif //KURSOVAYA_QUADRATICFORM_H
