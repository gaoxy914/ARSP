#include "point.h"

Point::Point() {
    m_dim = 0;
    m_coord = nullptr;
}

Point::Point(const int& dim, const double* coord) {
    m_dim = dim;
    m_coord = new double[m_dim];
    if (coord != nullptr) {
        memcpy(m_coord, coord, m_dim*sizeof(double));
    } else {
        memset(m_coord, 0, m_dim*sizeof(double));
    }
}

Point::Point(const Point& other) {
    m_dim = other.m_dim;
    m_coord = new double[m_dim];
    memcpy(m_coord, other.m_coord, m_dim*sizeof(double));
}

Point::~Point() {
    if (m_coord != nullptr) {
        delete[] m_coord;
        m_coord = nullptr;
    }
    m_dim = 0;
}

Point& Point::operator =(const Point& other) {
    if (&other != this) {
        if (m_dim != other.m_dim) {
            m_dim = other.m_dim;
            if (m_coord != nullptr) {
                delete[] m_coord;
            }
            m_coord = new double[m_dim];
        }
        memcpy(m_coord, other.m_coord, m_dim*sizeof(double));
    }
    return *this;
}

double Point::operator [](const int& index) const {
    assert(index >= 0 && index < m_dim);
    return m_coord[index];
}

bool operator ==(const Point& t, const Point& s) {
    if (t.m_dim != s.m_dim) return false;
    for (int i = 0; i < t.m_dim; ++ i) {
        if (t.m_coord[i] != s.m_coord[i]) {
            return false;
        }
    }
    return true;
}

ostream & operator <<(ostream& out, const Point& point) {
    out << "(";
    for (int i = 0; i < point.m_dim - 1; ++ i) {
        out << point[i] << ", ";
    }
    out << point[point.m_dim - 1] << ")";
    return out;
}