#include "hyperrect.h"

HyperRect::HyperRect() {
    m_dim = 0;
    m_lower = nullptr;
    m_upper = nullptr;
}

HyperRect::HyperRect(const int& dim, const double* lower, const double* upper) {
    m_dim = dim;
    m_lower = new double[m_dim];
    m_upper = new double[m_dim];
    if (lower != nullptr) {
        memcpy(m_lower, lower, m_dim*sizeof(double));
    } else {
        memset(m_lower, 0, m_dim*sizeof(double));
    }
    if (upper != nullptr) {
        memcpy(m_upper, upper, m_dim*sizeof(double));
    } else {
        memset(m_upper, 0, m_dim*sizeof(double));
    }
}

HyperRect::HyperRect(const Point& lower, const Point& upper) {
    m_dim = lower.m_dim;
    m_lower = new double[m_dim];
    m_upper = new double[m_dim];
    memcpy(m_lower, lower.m_coord, m_dim*sizeof(double));
    memcpy(m_upper, upper.m_coord, m_dim*sizeof(double));
}

HyperRect::HyperRect(const HyperRect& other) {
    m_dim = other.m_dim;
    m_lower = new double[m_dim];
    m_upper = new double[m_dim];
    memcpy(m_lower, other.m_lower, m_dim*sizeof(double));
    memcpy(m_upper, other.m_upper, m_dim*sizeof(double));
}

HyperRect::~HyperRect() {
    if (m_lower != nullptr) {
        delete[] m_lower;
        m_lower = nullptr;
    }
    if (m_upper != nullptr) {
        delete[] m_upper;
        m_upper = nullptr;
    }
    m_dim = 0;
}

HyperRect& HyperRect::operator =(const HyperRect& other) {
    if (&other != this) {
        if (m_dim != other.m_dim) {
            m_dim = other.m_dim;
            if (m_lower != nullptr) {
                delete[] m_lower;
            }
            m_lower = new double[m_dim];
            if (m_upper != nullptr) {
                delete[] m_upper;
            }
            m_upper = new double[m_dim];
        }
        memcpy(m_lower, other.m_lower, m_dim*sizeof(double));
        memcpy(m_upper, other.m_upper, m_dim*sizeof(double));
    }
    return *this;
}

void HyperRect::append(const Point& point) {
    assert(m_dim == point.m_dim);
    for (int i = 0; i < m_dim; ++ i) {
        m_lower[i] = min(m_lower[i], point[i]);
        m_upper[i] = max(m_upper[i], point[i]);
    }
}

void HyperRect::append(const HyperRect& rect) {
    assert(m_dim == rect.m_dim);
    for (int i = 0; i < m_dim; ++ i) {
        m_lower[i] = min(m_lower[i], rect.m_lower[i]);
        m_upper[i] = max(m_upper[i], rect.m_upper[i]);
    }
}

void HyperRect::append(const vector<double>& point) {
    assert(m_dim == point.size());
    for (int i = 0; i < m_dim; ++ i) {
        m_lower[i] = min(m_lower[i], point[i]);
        m_upper[i] = max(m_upper[i], point[i]);
    }
}

void HyperRect::get_vertices(vector<Point>& vertices) {
    int num = pow(2, m_dim);
    vertices.resize(num);
    for (int i = 0; i < num; ++ i) {
        vertices[i] = Point(m_dim + 1);
        for (int j = 0; j < m_dim; ++ j) {
            vertices[i].m_coord[j] = ((i>>j)&1) == 0 ? m_lower[j] : m_upper[j];
        }
        vertices[i].m_coord[m_dim] = 1;
    }
}

bool HyperRect::upper_edominate(const HyperRect& range) const {
    // cout << range << endl;
    bool unique = false;
    double min_value = -m_upper[m_dim - 1];
    for (int i = 0; i < m_dim - 1; ++ i) {
        if (m_upper[i] < 0) {
            unique = true;
            min_value -= m_upper[i]*range.m_lower[i];
        } else if (m_upper[i] > 0) {
            unique = true;
            min_value -= m_upper[i]*range.m_upper[i];
        }
    }
    return min_value >= 0 && unique;
}

bool HyperRect::lower_edominate(const HyperRect& range) const {
    bool unique = false;
    double min_value = -m_lower[m_dim - 1];
    for (int i = 0; i < m_dim - 1; ++ i) {
        if (m_lower[i] < 0) {
            unique = true;
            min_value -= m_lower[i]*range.m_lower[i];
        } else if (m_lower[i] > 0) {
            unique = true;
            min_value -= m_lower[i]*range.m_upper[i];
        }
    }
    return min_value >= 0 && unique;
}

HyperRect operator +(const HyperRect& rect1, const HyperRect& rect2) {
    assert(rect1.m_dim == rect2.m_dim);
    HyperRect rect(rect1);
    rect.append(rect2);
    return rect;
}

ostream & operator <<(ostream& out, const HyperRect& rect) {
    out << "[(";
    for (int i = 0; i < rect.m_dim - 1; ++ i) {
        out << rect.m_lower[i] << ", ";
    }
    out << rect.m_lower[rect.m_dim - 1] << ")\t(";
    for (int i = 0; i < rect.m_dim - 1; ++ i) {
        out << rect.m_upper[i] << ", ";
    }
    out << rect.m_upper[rect.m_dim - 1] << ")]";
    return out;
}