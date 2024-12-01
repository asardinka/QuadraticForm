//
// Created by Admin on 18.10.2024.
//

#ifndef KURSOVAYA_VECTOR_H
#define KURSOVAYA_VECTOR_H
#include <iostream>
using namespace std;

class Vector {
private:
    double* data;
    int size;

public:
    Vector(int s);
    ~Vector();

    double& operator[](int i);
    const double& operator[](int i) const;

    double dot(const Vector& other) const;
    double norm() const;
    int get_size() const;

    double product() const;
    double sum() const ;
    double sumOfSquares() const ;

    friend ostream& operator<<(ostream& s, const Vector& V);
};
#endif //KURSOVAYA_VECTOR_H
