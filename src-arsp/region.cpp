#include "region.h"

HalfSpace::HalfSpace() {
    m_dim = 0;
    m_coef = nullptr;
    m_side = true;
}

HalfSpace::HalfSpace(const int& dim, const double* coef, const bool& side) {
    m_dim = dim;
    m_side = side;
    m_coef = new double[m_dim + 1];
    if (coef != nullptr) {
        memcpy(m_coef, coef, (m_dim + 1)*sizeof(double));
    } else {
        memset(m_coef, 0, (m_dim + 1)*sizeof(double));
    }
}

HalfSpace::HalfSpace(const HalfSpace& other) {
    m_dim = other.m_dim;
    m_side = other.m_side;
    m_coef = new double[m_dim + 1];
    memcpy(m_coef, other.m_coef, (m_dim + 1)*sizeof(double));
}

HalfSpace::~HalfSpace() {
    if (m_coef != nullptr) {
        delete[] m_coef;
        m_coef = nullptr;
    }
    m_dim = 0;
    m_side = false;
}

HalfSpace& HalfSpace::operator =(const HalfSpace& other) {
    if (&other != this) {
        if (m_dim != other.m_dim) {
            m_dim = other.m_dim;
            if (m_coef != nullptr) {
                delete[] m_coef;
            }
            m_coef = new double[m_dim + 1];
        }
        memcpy(m_coef, other.m_coef, (m_dim + 1)*sizeof(double));
        m_side = other.m_side;
    }
    return *this;
}

ostream & operator <<(ostream& out, const HalfSpace& h) {
    for (int i = 0; i < h.m_dim - 1; ++ i) {
        out << h.m_coef[i] << " x w[" << i << "] + ";
    }
    out << h.m_coef[h.m_dim - 1] << " x w[" << h.m_dim - 1 << "]";
    out << (h.m_side ? " >= " : " <= ") << h.m_coef[h.m_dim];
    return out;
}

Region::Region(const int& dim, const int& c_cnt) : _inner(dim) {
    _lp = nullptr;
    _dim = dim;
    _c = c_cnt;
    _boundplane.resize(_c);
    for (int i = 0; i < _c; ++ i) {
        _boundplane[i] = HalfSpace(dim - 1); // w[d - 1] is replaced by 1 - w[0] ... - w[d - 2]
    }
    // _qpath_prefix = to_string(dim) + "_" + to_string(_c);
}

Region::~Region() {
    if(_lp != nullptr) {
        glp_delete_prob(_lp);
    }
}

double Region::_roundoff(const double& value, unsigned char prec) {
    double pow_10 = pow(10.0, (double)prec);
    return round(value*pow_10)/pow_10;
}

double Region::_rand_uniform(const double& a, const double& b) {
    random_device rd;
    default_random_engine generator(rd());
    uniform_real_distribution<double> distribution(a, b);
    return _roundoff(distribution(generator), PREC);
}

double Region::_rand_normal(const double& med, const double& var) {
    random_device rd;
    default_random_engine generator(rd());
    normal_distribution<double> distribution(med, var);
    return _roundoff(distribution(generator), PREC);
}

/* bool Region::compute_inner(const int& dim, double *inner) {
    // LP for computing inner point
    glp_set_obj_dir(m_lp, GLP_MIN);
    for (int i = 0; i < dim - 1; ++ i) {
        glp_set_obj_coef(m_lp, i + 1, 1);
    }
    glp_simplex(m_lp, nullptr);
    if (glp_get_status(m_lp) == GLP_OPT) {
        inner[dim - 1] = 1;
        for (int i = 0; i < dim - 1; ++ i) {
            inner[i] = glp_get_col_prim(m_lp, i + 1);
            inner[dim - 1] -= inner[i];
        }
        return true;
    } else {
        return false;
    }
#ifdef _DEBUG_
    assert(dim == m_inner.m_dim);
    assert(m_inner.m_coord != nullptr);
#endif
    for (int i = 0; i < dim; ++ i) {
        inner[i] = m_inner[i];
    }
    return true;
} */

Point& Region::get_inner() {
    return _inner;
}

Point& Region::get_vertex(const int& index) {
    return _vertex[index];
}

glp_prob* Region::get_lp() {
    return _lp;
}

vector<Point>& Region::get_vertex() {
    return _vertex;
}

int Region::get_vertex_size() {
    return _vertex.size();
}

/* void Region::gen_query() {
    srand((unsigned)time(nullptr));
    srand48((unsigned)time(nullptr));
    // generate a inner point
    double sum = 0;
    for (int i = 0; i < _dim; ++ i) {
        _inner.m_coord[i] = _rand_uniform(0, 1);
        sum += _inner[i];
    }
    for (int i = 0; i < _dim; ++ i) {
        _inner.m_coord[i] = _inner[i]/sum;
    }
#ifdef _DEBUG_
    cout << "inner point: " << _inner << endl;
#endif
    double p[_dim], alpha; // p = (1 - alpha)*d + alpha*inner
    double coef[_dim], beta;
    // generate constraints
    for (int i = 0; i < _c; ++ i) {
        alpha = _rand_uniform(0, 1);
        int d = rand()%_dim;
        for (int j = 0; j < _dim; ++ j) {
            p[j] = alpha*_inner[j];
            if (j == d) p[j] += (1 - alpha);
        }
        // a1(w1-p1) + ... + ad(1-w1-...-wd-1-pd)
        coef[_dim - 1] = _rand_uniform(0, 1);
        beta = coef[_dim - 1]*(_inner[_dim - 1] - p[_dim - 1]);
        _boundplane[i].m_coef[_dim - 1] = coef[_dim - 1]*(1 - p[_dim - 1]);
        for (int j = 0; j < _dim - 1; ++ j) {
            coef[j] = _rand_uniform(0, 1);
            beta += coef[j]*(_inner[j] - p[j]);
            _boundplane[i].m_coef[j] = (coef[j] - coef[_dim - 1]);
            _boundplane[i].m_coef[_dim - 1] -= coef[j]*p[j];
        }
        _boundplane[i].m_coef[_dim - 1] *= -1;
        _boundplane[i].m_side = beta >= 0;
#ifdef _DEBUG_
        cout << "generate the " << i << "-th constraints: \n";
        cout << "point = (";
        for (int j = 0; j < _dim - 1; ++ j) {
            cout << p[j] << ", ";
        }
        cout << p[_dim - 1] << ")\t";
        cout << "vector = (";
        for (int j = 0; j < _dim - 1; ++ j) {
            cout << coef[j] << ", ";
        }
        cout << coef[_dim - 1] << ")\t beta =" << beta << endl;
        cout << "after simplification: " << _boundplane[i] << endl;
#endif
    }
    write_query();
} */

void Region::gen_query() {
    srand((unsigned)time(nullptr));
    srand48((unsigned)time(nullptr));
    // generate a inner point
    double sum = 0;
    for (int i = 0; i < _dim; ++ i) {
        _inner.m_coord[i] = _rand_uniform(0, 1);
        sum += _inner[i];
    }
    for (int i = 0; i < _dim; ++ i) {
        _inner.m_coord[i] = _inner[i]/sum;
    }
    double t[_dim], s[_dim], side;
    for (int i = 0; i < _c; ++ i) {
        // generate two tuples
        side = 0;
        for (int j = 0; j < _dim; ++ j) {
            t[j] = drand48();
            s[j] = drand48();
            _boundplane[i].m_coef[j] = t[j] - s[j];
            side += _boundplane[i].m_coef[j]*_inner[j];
        }
        for (int j = 0; j < _dim - 1; ++ j) {
            _boundplane[i].m_coef[j] -= _boundplane[i].m_coef[_dim - 1];
        }
        _boundplane[i].m_coef[_dim - 1] *= -1;
        _boundplane[i].m_side = (side >= 0);
    }
    write_query();  
}

void Region::gen_weak_rankings() {
    assert(_c == _dim - 1);
    double sum = 0;
    for (int i = 0; i < _dim; ++ i) {
        _inner.m_coord[i] = _dim + 1 - i;
        sum += _inner[i];
    }
    for (int i = 0; i < _dim; ++ i) {
        _inner.m_coord[i] = _inner[i]/sum;
    }
    for (int i = 0; i < _c - 1; ++ i) {
        for (int j = 0; j < _dim; ++ j) {
            if (j == i) _boundplane[i].m_coef[j] = 1;
            if (j == i + 1) _boundplane[i].m_coef[j] = -1;
        }
        _boundplane[i].m_side = true;
    }
    for (int j = 0; j < _dim; ++ j) {
        if (j == _dim - 2) _boundplane[_c - 1].m_coef[j] = 2;
        else _boundplane[_c - 1].m_coef[j] = 1;
    }
    _boundplane[_c - 1].m_side = true;

    ofstream file((string("query/") + to_string(_dim) + string(".query")).c_str(), ios::out);
    for (int i = 0; i < _dim; ++ i) {
        file << _inner[i] << " ";
    }
    for (int i = 0; i < _c; ++ i) {
        for (int j = 0; j < _dim; ++ j) {
            file << _boundplane[i].m_coef[j] << " ";
        }
        file << _boundplane[i].m_side << " ";
    }
    file.close();
}

void Region::load_query() {
    ifstream file((string("query/") + to_string(_dim) + "_"\
     + to_string(_c) + string(".query")).c_str(), ios::in);
    for (int i = 0; i < _dim; ++ i) {
        file >> _inner.m_coord[i];
    }
    for (int i = 0; i < _c; ++ i) {
        for (int j = 0; j < _dim; ++ j) {
            file >> _boundplane[i].m_coef[j];
        }
        file >> _boundplane[i].m_side;
    }
    file.close();
#ifdef _DEBUG_
    printf("Query Parameters: dim = %d, c = %d.\n", _dim, _c);
    cout << _inner << endl;
    for (int i = 0; i < _c; ++ i) {
        cout << _boundplane[i] << endl;
    }
#endif
}

void Region::load_weak_rankings() {
    ifstream file((string("query/") + to_string(_dim) + string(".query")).c_str(), ios::in);
    for (int i = 0; i < _dim; ++ i) {
        file >> _inner.m_coord[i];
    }
    for (int i = 0; i < _c; ++ i) {
        for (int j = 0; j < _dim; ++ j) {
            file >> _boundplane[i].m_coef[j];
        }
        file >> _boundplane[i].m_side;
    }
    file.close();
#ifdef _DEBUG_
    printf("Query Parameters: dim = %d, c = %d.\n", _dim, _c);
    cout << _inner << endl;
    for (int i = 0; i < _c; ++ i) {
        cout << _boundplane[i] << endl;
    }
#endif
}

void Region::load_nba_query() {
    ifstream file("query/nba.query", ios::in);
    for (int i = 0; i < _dim; ++ i) {
        file >> _inner.m_coord[i];
    }
    for (int i = 0; i < _c; ++ i) {
        for (int j = 0; j < _dim; ++ j) {
            file >> _boundplane[i].m_coef[j];
        }
        file >> _boundplane[i].m_side;
    }
    file.close();
    printf("Query Parameters: dim = %d, c = %d.\n", _dim, _c);
    cout << _inner << endl;
    for (int i = 0; i < _c; ++ i) {
        cout << _boundplane[i] << endl;
    }
}

void Region::write_query() {
    ofstream file((string("query/") + to_string(_dim) + "_"\
     + to_string(_c) + string(".query")).c_str(), ios::out);
    for (int i = 0; i < _dim; ++ i) {
        file << _inner[i] << " ";
    }
    for (int i = 0; i < _c; ++ i) {
        for (int j = 0; j < _dim; ++ j) {
            file << _boundplane[i].m_coef[j] << " ";
        }
        file << _boundplane[i].m_side << " ";
    }
    file.close();
#ifdef _DEBUG_
    printf("Query Parameters: dim = %d, c = %d.\n", _dim, _c);
    cout << _inner << endl;
    for (int i = 0; i < _c; ++ i) {
        cout << _boundplane[i] << endl;
    }
#endif
}

void Region::add_plane() {
    srand((unsigned)time(nullptr));
    srand48((unsigned)time(nullptr));
    double t[_dim], s[_dim], side;
    // generate two tuples
    side = 0;
    HalfSpace plane(_dim - 1);
    for (int j = 0; j < _dim; ++ j) {
        t[j] = drand48();
        s[j] = drand48();
        plane.m_coef[j] = t[j] - s[j];
        side += plane.m_coef[j]*_inner[j];
    }
    for (int j = 0; j < _dim - 1; ++ j) {
        plane.m_coef[j] -= plane.m_coef[_dim - 1];
    }
    plane.m_coef[_dim - 1] *= -1;
    plane.m_side = (side >= 0);
    _boundplane.push_back(plane);
    _c ++;
    write_query();
}

void Region::build_lp() {
#ifdef _DEBUG_
    assert(!_boundplane.empty());
#endif
    _lp = glp_create_prob();
    glp_term_out(GLP_OFF);
    glp_add_cols(_lp, _dim - 1);
    for (int i = 0; i < _dim - 1; ++ i) {
        glp_set_col_bnds(_lp, i + 1, GLP_DB, 0.0, 1.0);
    }
    int c = _boundplane.size();
    glp_add_rows(_lp, c + 1);
    int size = (_dim - 1)*(c + 1);
    for (int i = 0; i < _boundplane.size(); ++ i) {
        if (_boundplane[i].m_side) { // Ax >= b
            glp_set_row_bnds(_lp, i + 1, GLP_LO, _boundplane[i].m_coef[_dim - 1], 0.0);
        } else {
            glp_set_row_bnds(_lp, i + 1, GLP_UP, 0.0, _boundplane[i].m_coef[_dim - 1]);
        }
    }
    glp_set_row_bnds(_lp, c + 1, GLP_DB, 0.0, 1.0);
    int ia[size + 1];
    int ja[size + 1];
    double ar[size + 1];
    for (int i = 0; i < c + 1; ++ i) {
        for (int j = 0; j < _dim - 1; ++ j) {
            int k = i*(_dim - 1) + j + 1;
            ia[k] = i + 1;
            ja[k] = j + 1;
            if (i == c) { ar[k] = 1; }
            else { ar[k] = _boundplane[i].m_coef[j]; }
        }
    }
    glp_load_matrix(_lp, size, ia, ja, ar);
}

void Region::compute_vertex() {
#ifdef _DEBUG_
    assert(!_boundplane.empty());
#endif
    if (_dim > 2) {
        // Ax + b <= 0
        int n = _c + 2*(_dim - 1) + 1;
        int size = _dim*n;
        double normals[size];
        int index = 0;
        for (int i = 0; i < _dim; ++ i) { // x1+...+xd-1 <= 1
            if (i == _dim - 1) normals[index ++] = -1;
            else if (i < _dim - 1) normals[index ++] = 1;
        }
        for (int i = 0; i < _dim - 1; ++ i) {
            for (int j = 0; j < _dim; ++ j) { // xi <= 1
                if (i == j) normals[index ++] = 1;
                else normals[index ++] = 0;
            }
            normals[index - 1] = -1;
            for (int j = 0; j < _dim; ++ j) { // xi >= 0
                if (i == j) normals[index ++] = -1;
                else normals[index ++] = 0;
            }
        }
        for (int i = 0; i < _c; ++ i) {
            if (_boundplane[i].m_side) {
                for (int j = 0; j < _dim; ++ j) {
                    if (j == _dim - 1) normals[index ++] = _boundplane[i].m_coef[j];
                    else if (j < _dim - 1) normals[index ++] = -_boundplane[i].m_coef[j];
                }
            } else {
                for (int j = 0; j < _dim; ++ j) {
                    if (j == _dim - 1) normals[index ++] = -_boundplane[i].m_coef[j];
                    else if (j < _dim - 1) normals[index ++] = _boundplane[i].m_coef[j];
                }
            }
        }
        Coordinates feasible;
        for (int i = 0; i < _dim - 1; ++ i) feasible << _inner[i];
        // cout << feasible << endl;
        Qhull q;
        q.setFeasiblePoint(feasible);
        q.runQhull("", _dim, n, normals, "H");
#ifdef _DEBUG_
        cout << q.facetCount() << endl;
#endif
        _vertex.resize(q.facetCount());
        QhullFacetListIterator it(q.facetList());
        index = 0;
        while (it.hasNext()) {
            QhullHyperplane plane = it.next().hyperplane();
            _vertex[index] = Point(_dim);
            _vertex[index].m_coord[_dim - 1] = 1;
            for (int i = 0; i < _dim - 1; ++ i) {
                _vertex[index].m_coord[i] = -plane[i]/plane.offset() + _inner[i];
                if (_vertex[index].m_coord[i] < 0) _vertex[index].m_coord[i] = 0;
                _vertex[index].m_coord[_dim - 1] -= _vertex[index][i];
            }
            if (_vertex[index].m_coord[_dim - 1] < 0) _vertex[index].m_coord[_dim - 1] = 0;
            ++ index;
        }
    } else if (_dim == 2) {
        double left = 0, right = 1;
        for (int i = 0; i < _boundplane.size(); ++ i) {
            if (_boundplane[i].m_coef[0] == 0) continue;
            double beta = _boundplane[i].m_coef[1]/_boundplane[i].m_coef[0];
            if (_boundplane[i].m_side) {
                if (_boundplane[i].m_coef[0] > 0) {
                    left = max(left, beta);
                } else if (_boundplane[i].m_coef[0] < 0) {
                    right = min(right, beta);
                }
            } else {
                if (_boundplane[i].m_coef[0] > 0) {
                    right = min(right, beta);
                } else if (_boundplane[i].m_coef[0] < 0) {
                    left = max(left, beta);
                }
            }
#ifdef _DEBUG_
        assert(left <= right);
#endif
        _vertex.resize(2);
        _vertex[0] = Point(2);
        _vertex[0].m_coord[0] = left;
        _vertex[0].m_coord[1] = 1 - left;
        _vertex[1] = Point(2);
        _vertex[1].m_coord[0] = right;
        _vertex[1].m_coord[1] = 1 - right;
        }
    }
#ifdef _DEBUG_
    cout << _vertex.size() << endl;
    for (int i = 0; i < _vertex.size(); ++ i) {
        cout << _vertex[i] << endl;
    }
#endif
}

void Region::print() {
    printf("Query Parameters: dim = %d, c = %d.\n", _dim, _c);
    cout << _inner << endl;
    for (int i = 0; i < _c; ++ i) {
        cout << _boundplane[i] << endl;
    }
    cout << _vertex.size() << endl;
    for (int i = 0; i < _vertex.size(); ++ i) {
        cout << _vertex[i] << endl;
    }
}