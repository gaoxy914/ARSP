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

double HyperRect::spherical_volume() const {
    double sum_of_sqaures = 0;
    double radius;
    for (int i = 0; i < m_dim; ++ i) {
        double half_side_length = (m_upper[i] - m_lower[i])*0.5;
        sum_of_sqaures += half_side_length*half_side_length;
    }
    radius = sqrt(sum_of_sqaures);
    if (m_dim == 3) {
        return (radius*radius*radius*UNIT_SPHERE_VOLUMES[m_dim]);
    } else if (m_dim == 2) {
        return (radius*radius*UNIT_SPHERE_VOLUMES[m_dim]);
    } else {
        return (pow(radius, m_dim)*UNIT_SPHERE_VOLUMES[m_dim]);
    }
}

double HyperRect::lower_score(const Point& weight) const {
    double s = weight[m_dim - 1]*m_lower[m_dim - 1];
    for (int i = 0; i < m_dim - 1; ++ i) {
        s += weight[i]*m_lower[i];
    }
    return s;
}

inline double _dot(const int& dim, const double* t, const double* s) {
    double sum = 0;
    for (int i = 0; i < dim; ++ i) {
        sum += t[i]*s[i];
    }
    return sum;
}

bool HyperRect::upper_fdominate_V(const Point& point, const vector<Point>& vertex) const {
    bool unique = false;
    double score_t = 0, score_s = 0;
    for (int i = 0; i < vertex.size(); ++ i) {
        score_t = _dot(m_dim, m_upper, vertex[i].m_coord);
        score_s = _dot(m_dim, point.m_coord, vertex[i].m_coord);
        if (score_t > score_s) {
            return false;
        } else if (score_t < score_s) {
            unique = true;
        }
    }
    return unique;
}

bool HyperRect::lower_fdominate_V(const Point& point, const vector<Point>& vertex) const {
    bool unique = false;
    double score_t = 0, score_s = 0;
    for (int i = 0; i < vertex.size(); ++ i) {
        score_t = _dot(m_dim, m_lower, vertex[i].m_coord);
        score_s = _dot(m_dim, point.m_coord, vertex[i].m_coord);
        if (score_t > score_s) {
            return false;
        } else if (score_t < score_s) {
            unique = true;
        }
    }
    return unique;
}

bool HyperRect::upper_edominate(const HyperRect& range) const {
    // cout << range << endl;
    double min_value = -m_upper[m_dim - 1];
    for (int i = 0; i < m_dim - 1; ++ i) {
        if (m_upper[i] < 0) {
            min_value -= m_upper[i]*range.m_lower[i];
        } else {
            min_value -= m_upper[i]*range.m_upper[i];
        }
    }
    return min_value >= 0;
}

bool HyperRect::lower_edominate(const HyperRect& range) const {
    double min_value = -m_lower[m_dim - 1];
    for (int i = 0; i < m_dim - 1; ++ i) {
        if (m_lower[i] < 0) {
            min_value -= m_lower[i]*range.m_lower[i];
        } else {
            min_value -= m_lower[i]*range.m_upper[i];
        }
    }
    return min_value >= 0;
}
    

void HyperRect::z_split(vector<HyperRect>& rect, Point& mid_point) {
#ifdef _DEBUG_   
    assert(mid_point.m_dim == m_dim);
#endif
    int cnt = pow(2, m_dim);
    rect.reserve(cnt);
    for (int i = 0; i < m_dim; ++ i) {
        mid_point.m_coord[i] = (m_lower[i] + m_upper[i])/2;
    }
    for (int i = 0; i < cnt; ++ i) {
        HyperRect sub_rect = HyperRect(m_dim);
        for (int j = 0; j < m_dim; ++ j) {
            sub_rect.m_lower[j] = (i>>j)&1 ? mid_point[j] : m_lower[j];
            sub_rect.m_upper[j] = (i>>j)&1 ? m_upper[j] : mid_point[j];
        }
        rect.push_back(sub_rect);
    }
}

HyperRect HyperRect::get_subrect(const int& k) {
    HyperRect sub_rect(m_dim);
    for (int i = 0; i < m_dim; ++ i) {
        sub_rect.m_lower[i] = (k>>i)&1 ? (m_lower[i] + m_upper[i])/2 : m_lower[i];
        sub_rect.m_upper[i] = (k>>i)&1 ? m_upper[i] : (m_lower[i] + m_upper[i])/2;
    }
    return sub_rect;
}

Point HyperRect::get_mid() {
    Point mid_point(m_dim);
    for (int i = 0; i < m_dim; ++ i) {
        mid_point.m_coord[i] = (m_lower[i] + m_upper[i])/2;
    }
    return mid_point;
}

// plane: w[1]x[1] + ... + w[d-1]x[d-1] = w[d]
bool HyperRect::intersect(const vector<double>& plane) const {
    assert(m_dim == plane.size() - 1);
    double min_val = 0, max_val = 0;
    for (int i = 0; i < m_dim; ++ i) {
        if (plane[i] < 0) {
            min_val += m_upper[i]*plane[i];
            max_val += m_lower[i]*plane[i];
        } else {
            min_val += m_lower[i]*plane[i];
            max_val += m_upper[i]*plane[i];
        }
    }
    return min_val <= plane[m_dim] && max_val >= plane[m_dim];
}

bool HyperRect::contain(const HyperRect& other) const {
    for (int i = 0; i < m_dim; ++ i) {
        if (m_lower[i] > other.m_lower[i] || m_upper[i] < other.m_upper[i]) {
            return false;
        }
    }
    return true;
}

bool HyperRect::intersect(const HyperRect& other) const {
    for (int i = 0; i < m_dim; ++ i) {
        if (m_lower[i] > other.m_upper[i] || m_upper[i] < other.m_lower[i]) {
            return false;
        }
    }
    return true;
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