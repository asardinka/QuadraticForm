#include "Vector.h"
#include  "math.h"
Vector::Vector(int s) : size(s) {
    data = new double[size];
}

Vector::~Vector() {
    delete[] data;
}

double& Vector::operator[](int i) {
    return data[i];
}

const double& Vector::operator[](int i) const {
    return data[i];
}

double Vector::dot(const Vector& other) const {
    double result = 0;
    for (int i = 0; i < size; ++i) {
        result += (data[i] * other[i]);
    }
    return result;
}

double Vector::norm() const {
    return double(sqrt(dot(*this)));
}

int Vector::get_size() const {
    return size;
}

ostream& operator<<(std::ostream& s, const Vector& V) {
    for (int i = 0; i < V.size; i++) {
        s << V.data[i]  << std::endl;
    }
    return s;
}

double Vector::product() const {
    double result = 1.0;
    for (int i = 0; i < size; ++i) {
        result *= data[i];
    }
    return result;
}

double Vector::sum() const {
    double result = 0.0;
    for (int i = 0; i < size; ++i) {
        result += data[i];
    }
    return result;
}

double Vector::sumOfSquares() const {
    double result = 0.0;
    for (int i = 0; i < size; ++i) {
        result += data[i] * data[i];
    }
    return result;
}