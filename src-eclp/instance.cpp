#include "instance.h"

Instance::Instance() {
    m_id = -1;
    m_obj_id = -1;
}

Instance::Instance(const int& dim, const double* coord, const int& id, const int& obj_id) : m_point(dim, coord) {
    m_id = id;
    m_obj_id = obj_id;
}

Instance::Instance(const Instance& other) : m_point(other.m_point) {
    m_id = other.m_id;
    m_obj_id = other.m_obj_id;
}

Instance::~Instance() {
    m_id = -1;
    m_obj_id = -1;
}

double Instance::score(const Point& weight) const {
    int dim = weight.m_dim;
    double s = weight[dim - 1]*m_point[dim - 1];
    for (int i = 0; i < dim - 1; ++ i) {
        s += weight[i]*m_point[i];
    }
    return s;
}

bool Instance::dominate(const HyperRect& rect) const {
    bool unique = false;
    for (int i = 0; i < m_point.m_dim; ++ i) {
        if (m_point[i] > rect.m_lower[i]) {
            return false;
        } else if (m_point[i] < rect.m_lower[i]) {
            unique = true;
        }
    }
    return unique;
}

bool Instance::dominate(const Instance& other) const {
    bool unique = false;
    for (int i = 0; i < m_point.m_dim; ++ i) {
        if (m_point[i] > other.m_point[i]) {
            return false;
        } else if (m_point[i] < other.m_point[i]) {
            unique = true;
        }
    }
    return unique;
}

bool Instance::cross(const HyperRect& rect) const {
    bool unique = false;
    for (int i = 0; i < m_point.m_dim; ++ i) {
        if (m_point[i] > rect.m_upper[i]) {
            return false;
        } else if (m_point[i] < rect.m_upper[i]) {
            unique = true;
        }
    }
    return unique;
}

Instance& Instance::operator =(const Instance& other) {
    if (&other != this) {
        m_point = other.m_point;
        m_id = other.m_id;
        m_obj_id = other.m_obj_id;
    }
    return *this;
}

ostream & operator <<(ostream& out, const Instance& t) {
    out << "id = " << t.m_id << ", obj_id = " << t.m_obj_id << ", ";
    out << t.m_point;
    return out;
}