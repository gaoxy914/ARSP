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

Instance& Instance::operator =(const Instance& other) {
    if (&other != this) {
        m_point = other.m_point;
        m_id = other.m_id;
        m_obj_id = other.m_obj_id;
    }
    return *this;
}

double Instance::score(const int& dim, const double* weight) const {
#ifdef _DEBUG_
    assert(dim == m_point.m_dim);
#endif
    double s = weight[dim - 1]*m_point[dim - 1];
    for (int i = 0; i < dim - 1; ++ i) {
        s += weight[i]*m_point[i];
    }
    return s;
}

double Instance::score(const Point& weight) const {
    return score(weight.m_dim, weight.m_coord);
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

bool Instance::fdominate_LP(const Instance& other, glp_prob *lp) const {
#ifdef _DEBUG_
    assert(lp);
    if (other.m_id == 2) cout << m_id << '\t' << other.m_id << endl;
#endif
    glp_set_obj_dir(lp, GLP_MIN);
    int dim = m_point.m_dim;
    double coef, delta = other.m_point[dim - 1] - m_point[dim - 1];
    for (int i = 0; i < dim - 1; ++ i) {
        coef = other.m_point[i] - m_point[i] - delta;
        glp_set_obj_coef(lp, i + 1, coef);
    }
    glp_simplex(lp, nullptr);
#ifdef _DEBUG_
    assert(glp_get_status(lp) == GLP_OPT);
#endif
    if (glp_get_status(lp) == GLP_OPT) {
        return glp_get_obj_val(lp) + delta >= 0;
    } else {
        return false;
    }
}

bool Instance::fdominate_LP(const Instance& other, glp_prob *lp, Point& pivot) const {
    glp_set_obj_dir(lp, GLP_MIN);
    int dim = m_point.m_dim;
    double coef, delta = other.m_point[dim - 1] - m_point[dim - 1];
    for (int i = 0; i < dim - 1; ++ i) {
        coef = other.m_point[i] - m_point[i] - delta;
        glp_set_obj_coef(lp, i + 1, coef);
    }
    glp_simplex(lp, nullptr);
#ifdef _DEBUG_
    assert(glp_get_status(lp) == GLP_OPT);
#endif
    if (glp_get_status(lp) == GLP_OPT) {
        return glp_get_obj_val(lp) + delta >= 0;
    } else {
        return false;
    }
}

bool Instance::fdominate_V(const Instance& other, const vector<Point>& vertex) const {
    bool unique = false;
    double score_t = 0, score_s = 0;
    for (int i = 0; i < vertex.size(); ++ i) {
        score_t = score(vertex[i]);
        score_s = other.score(vertex[i]);
        if (score_t > score_s) {
            return false;
        } else if (score_t < score_s) {
            unique = true;
        }
    }
    return unique;
}

bool Instance::fdominate_V(const HyperRect& rect, const vector<Point>& vertex) const {
    bool unique =false;
    double score_t = 0, score_s = 0;
    for (int i = 0; i < vertex.size(); ++ i) {
        score_t = score(vertex[i]);
        score_s = rect.lower_score(vertex[i]);
        if (score_t > score_s) {
            return false;
        } else if (score_t < score_s) {
            unique = true;
        }
    }
    return unique;
}

bool Instance::edominate(const Instance& other, const HyperRect& rect) const {
    bool unique = false;
    int dim = m_point.m_dim;
    double min_value = other.m_point[dim - 1] - m_point[dim - 1];
    for (int i = 0; i < dim - 1; ++ i) {
        if (other.m_point[i] > m_point[i]) {
            unique = true;
            min_value += (other.m_point[i] - m_point[i])*rect.m_lower[i];
        } else if (other.m_point[i] < m_point[i]) {
            unique = true;
            min_value += (other.m_point[i] - m_point[i])*rect.m_upper[i];
        }
    }
    return min_value >= 0 && unique;
}

ostream & operator <<(ostream& out, const Instance& t) {
    out << "id = " << t.m_id << ", obj_id = " << t.m_obj_id << ", ";
    out << t.m_point;
    return out;
}