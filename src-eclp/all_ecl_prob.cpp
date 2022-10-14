#include "all_ecl_prob.h"

double Dataset::_roundoff(const double& value, unsigned char prec) {
    double pow_10 = pow(10.0, (double)prec);
    return round(value*pow_10)/pow_10;
}

double Dataset::_rand_uniform(const double& a, const double& b) {
    return drand48()*(b - a) + a;
}

double Dataset::_rand_normal(const double& med, const double& var) {
    double x, y, s;
    do {
        x = drand48();
        if (rand()%2) x = -x;
        y = drand48();
        if (rand()%2) y = -y;
        s = x*x + y*y;
    } while (s >= 1);
    double num = x*sqrt(-2.*log(s)/s);
    num = num*var + med;
    return num;
}

void Dataset::_gen_inde_data() {
    double center[_dim], lower[_dim], upper[_dim];
    int start = 0;
    double length = 0;
    bool phi_flag = false;
    for (int i = 0; i < _m; ++ i) {
        for (int j = 0; j < _dim; ++ j)
            center[j] = _rand_uniform(0, 1);
        for (int j = 0; j < _dim; ++ j) {
            do { 
                length = _rand_normal(_l/2.0, _l/8.0);
            } while (length > _l || length <= 0);
            lower[j] = center[j] - length/2 > 0 ? center[j] - length/2 : 0;
            upper[j] = center[j] + length/2 < 1 ? center[j] + length/2 : 1;
        }
        if (_rand_uniform(0, 1) < _phi_prob) {
            phi_flag = true;
        }
        for (int j = start; j < start + _obj_cnt[i]; ++ j) {
            _data[j].m_obj_id = i;
            _data[j].m_id = j;
            _data[j].m_point = Point(_dim);
            for (int k = 0; k < _dim; ++ k) {
                _data[j].m_point.m_coord[k] = _rand_uniform(lower[k], upper[k]);
            }
        }
        start += _obj_cnt[i];
        if (phi_flag) { // add empty instance to object
            for (int j = 0; j < _dim; ++ j) {
                _data[start - 1].m_point.m_coord[j] = 1;
            }
        }
        phi_flag = false;
    }
}

void Dataset::_gen_anti_data() {
    double c = 0;
    double center[_dim], lower[_dim], upper[_dim];
    int start = 0;
    double length = 0;
    bool phi_flag = false;
    for (int i = 0; i < _m; ++ i) {
        do {
            c = _rand_normal(0.5, 0.05);
        } while (c >= 1 || c <= 0);
        do {
            center[_dim - 1] = _dim*c;
            for (int j = 0; j < _dim - 1; ++ j) {
                center[j] = _rand_uniform(0, 1);
                center[_dim - 1] -= center[j];
            }
        } while (center[_dim - 1] < 0 || center[_dim - 1] > 1);
        for (int j = 0; j < _dim; ++ j) {
            do { 
                length = _rand_normal(_l/2.0, _l/8.0);
            } while (length > _l || length <= 0);
            lower[j] = center[j] - length/2 > 0 ? center[j] - length/2 : 0;
            upper[j] = center[j] + length/2 < 1 ? center[j] + length/2 : 1;
        }
        if (_rand_uniform(0, 1) < _phi_prob) {
            phi_flag = true;
        }
        for (int j = start; j < start + _obj_cnt[i]; ++ j) {
            _data[j].m_obj_id = i;
            _data[j].m_id = j;
            _data[j].m_point = Point(_dim);
            for (int k = 0; k < _dim; ++ k) {
                _data[j].m_point.m_coord[k] = _rand_uniform(lower[k], upper[k]);
            }
        }
        start += _obj_cnt[i];
        if (phi_flag) { // add empty instance to object
            for (int j = 0; j < _dim; ++ j) {
                _data[start - 1].m_point.m_coord[j] = 1;
            }
        }
        phi_flag = false;
    } 
}

void Dataset::_gen_corr_data() {
    double c = 0;
    double center[_dim], lower[_dim], upper[_dim];
    int start = 0;
    double length = 0;
    bool phi_flag = false;
    for (int i = 0; i < _m; ++ i) {
        do {
            c = _rand_normal(0.5, 0.05);
        } while (c >= 1 || c <= 0);
        do {
            center[_dim - 1] = _dim*c;
            for (int j = 0; j < _dim - 1; ++ j) {
                center[j] = _rand_normal(c, 0.05);
                center[_dim - 1] -= center[j];
            }
        } while (center[_dim - 1] < 0 || center[_dim - 1] > 1);
        for (int j = 0; j < _dim; ++ j) {
            do {
                length = _rand_normal(_l/2.0, _l/8.0);
            } while (length > _l || length <= 0);
            lower[j] = center[j] - length/2 > 0 ? center[j] - length/2 : 0;
            upper[j] = center[j] + length/2 < 1 ? center[j] + length/2 : 1;
        }
        if (_rand_uniform(0, 1) < _phi_prob) {
            phi_flag = true;
        }
        for (int j = start; j < start + _obj_cnt[i]; ++ j) {
            _data[j].m_obj_id = i;
            _data[j].m_id = j;
            _data[j].m_point = Point(_dim);
            for (int k = 0; k < _dim; ++ k) {
                _data[j].m_point.m_coord[k] = _rand_uniform(lower[k], upper[k]);
            }
        }
        start += _obj_cnt[i];
        if (phi_flag) { // add empty instance to object
            for (int j = 0; j < _dim; ++ j) {
                _data[start - 1].m_point.m_coord[j] = 1;
            }
        }
        phi_flag = false;
    }
}

int Dataset::_filter_zero_skyprob() {
    unordered_map<int, double> skyprob;
    KDTree kdtree(_dim);
    kdtree.load_data(_n, _data);
    kdtree.mix_sky_prob(_m, _obj_cnt, skyprob);
    for (int i = 0; i < _n; ++ i) {
        if (skyprob.count(_data[i].m_id) == 0) {
            swap(_data[i], _data[_n - 1]);
            -- _n;
            -- i;
        }
    }
    assert(_n == skyprob.size());
    cout << "#instance after preprcessing: " << _n << endl;
    return _n;
}

Dataset::Dataset(const int& dim, const int& m, const int& cnt, const double& l, const double& phi_prob) {
    _dim = dim;
    _m = m;
    _cnt = cnt;
    _l = l;
    _phi_prob = phi_prob;

    _n = 0;
    _obj_cnt = new int[m];
    _data = nullptr;
    _dpath_prefix = to_string(_dim) + "_" + to_string(_m) + "_" + to_string(_cnt);
    stringstream stream;
    stream << fixed << setprecision(1) << _l << "_" << fixed << setprecision(2) << _phi_prob;
    _dpath_prefix += "_" + stream.str();

    srand((unsigned)time(nullptr));
    srand48((unsigned)time(nullptr));
}

Dataset:: ~Dataset() {
    if (_obj_cnt != nullptr) {
        delete[] _obj_cnt;
    }
    if (_data != nullptr) {
        delete[] _data;
    }
}

void Dataset::gen_data(const char *path, const int& type) {
    for (int i = 0; i < _m; ++ i) {
        _obj_cnt[i] = rand()%_cnt + 1;
        _n += _obj_cnt[i];
    }
    _data = new Instance[_n];
    switch (type) {
    case 1:
        _gen_inde_data();
        break;
    case 2:
        _gen_anti_data();
        break;
    case 3:
        _gen_corr_data();
        break;
    default:
        break;
    }
    write_data(path);
}

void Dataset::load_data(const char *path) {
    ifstream file((string(path) + _dpath_prefix + string(".cnt")).c_str(), ios::in);
    if (!file.is_open()) {
        printf("Fail in loading data.\n");
        exit(1);
    }
    for (int i = 0; i < _m; ++ i) {
        file >> _obj_cnt[i];
        _n += _obj_cnt[i];
    }
    file.close();
    file.open((string(path) + _dpath_prefix + string(".data")).c_str(), ios::in);
    _data = new Instance[_n];
    for (int i = 0; i < _n; ++ i) {
        file >> _data[i].m_obj_id;
        file >> _data[i].m_id;
        _data[i].m_point = Point(_dim);
        for (int j = 0; j < _dim; ++ j) {
            file >> _data[i].m_point.m_coord[j];
        }
    }
    file.close();
#ifdef _DEBUG_
    printf("Dataset Parameters: dim = %d, m = %d, cnt = %d, l = %.1f,\
     phi = %.1f, n = %d.\n", _dim, _m, _cnt, _l, _phi_prob, _n);
    for (int i = 0; i < _n; ++ i) {
        cout << _data[i] << endl;
    }
#endif
}

void Dataset::write_data(const char *path) {
    ofstream file((string(path) + _dpath_prefix + string(".cnt")).c_str(), ios::out);
    for (int i = 0; i < _m; ++ i) {
        file << _obj_cnt[i] << " ";
    }
    file.close();
    file.open((string(path) + _dpath_prefix + string(".data")).c_str(), ios::out);
    for (int i = 0; i < _n; ++ i) {
        file << _data[i].m_obj_id << " ";
        file << _data[i].m_id << " ";
        for (int j = 0; j < _dim; ++ j) {
            file << _data[i].m_point[j] << " ";
        }
    }
    file.close();
#ifdef _DEBUG_
    printf("Dataset Parameters: dim = %d, m = %d, cnt = %d, l = %.1f,\
     phi = %.1f, n = %d.\n", _dim, _m, _cnt, _l, _phi_prob, _n);
    for (int i = 0; i < _n; ++ i) {
        cout << _data[i] << endl;
    }
#endif
}

void Dataset::load_nba_data() {
    ifstream file("data/nba.cnt", ios::in);
    if (!file.is_open()) {
        printf("Fail in loading data.\n");
        exit(1);
    }
    for (int i = 0; i < _m; ++ i) {
        file >> _obj_cnt[i];
        _n += _obj_cnt[i];
    }
    file.close();
    cout << _n << endl;
    file.open("data/nba.data", ios::in);
    _data = new Instance[_n];
    int obj_id = -1;
    int cnt = 0;
    double temp;
    for (int i = 0; i < _n; ++ i) {
        if (i == cnt) {
            ++ obj_id;
            cnt += _obj_cnt[obj_id];
        }
        _data[i].m_obj_id = obj_id;
        _data[i].m_id = i;
        _data[i].m_point = Point(_dim);
        for (int j = 0; j < 5; ++ j) {
            if (j < _dim) file >> _data[i].m_point.m_coord[j];
            else file >> temp;
        }
    }
    file.close();
    /* for (int i = 0; i < _n; ++ i) {
        cout << _data[i] << endl;
    } */
}

#define mp make_pair

void Dataset::dual_ms_preprocess_2d() {
    _filter_zero_skyprob();

    rsky_prob.resize(_n);
    Theta.resize(_n);
    
    int sigma[_m];
    double beta, delta, x, y, theta;
    // test for point 0
    for (int i = 0; i < _n; ++ i) {
        // handle point i, compute all possible cases
        // cout << i << endl;
        memset(sigma, 0, _m*sizeof(int));
        beta = 1;
        vector<pair<double, int> > theta2id;
        for (int j = 0; j < _n; ++ j) {
            if (_data[i].m_obj_id != _data[j].m_obj_id) {
                if (_data[j].dominate(_data[i])) {
                    delta = _obj_cnt[_data[j].m_obj_id] - sigma[_data[j].m_obj_id];
                    beta *= (delta - 1)/delta;
                    ++ sigma[_data[j].m_obj_id];
                } else if (_data[i].dominate(_data[j])) {
                    continue;
                } else {
                    x = _data[j].m_point[0] - _data[i].m_point[0];
                    y = _data[j].m_point[1] - _data[i].m_point[1];
                    theta = _roundoff(atan2(y, x)*180/PI, PREC);
                    if (y < 0) theta += 360;
                    theta2id.push_back(mp(theta, _data[j].m_obj_id));
                    Theta[i].push_back(theta);
                }
            }
        }

        // sort according to theta
        sort(theta2id.begin(), theta2id.end(), [](const pair<double, int>& p1, const pair<double, int>& p2){
            return p1.first < p2.first;
        });
        sort(Theta[i].begin(), Theta[i].end());

        for (int s = 0; s < theta2id.size(); ++ s) {
            for (int e = s; e < theta2id.size(); ++ e) {
                int index = s;
                // compute prob
                for ( ; index <= e; ++ index) {
                    delta = _obj_cnt[theta2id[index].second] - sigma[theta2id[index].second];
                    if (delta == 1) {
                        break;
                    } else if (delta > 1) {
                        beta *= (delta - 1)/delta;
                        ++ sigma[theta2id[index].second];
                    }
                }
                if (index == e + 1) {
                    rsky_prob[i].insert(mp(mp(theta2id[s].first, theta2id[e].first), beta*1.0/_obj_cnt[_data[i].m_obj_id]));
                } else {
                    e = theta2id.size();
                }
                // recover sigma
                -- index;
                for ( ; index >= s; -- index) {
                    delta = _obj_cnt[theta2id[index].second] - sigma[theta2id[index].second];
                    beta *= (delta + 1)/delta;
                    -- sigma[theta2id[index].second];
                }
            }
        }
    }
}

void Dataset::dual_ms_2d(HyperRect& ratios, unordered_map<int, double>& result) {
    /* for (int i = 0; i < _n; ++ i) {
        cout << _data[i] << endl;
    } */
    double alpha = _roundoff(atan2(-ratios.m_lower[0], 1)*180/PI, PREC) + 180;
    double beta = _roundoff(atan2(-ratios.m_upper[0], 1)*180/PI, PREC) + 360;
    for (int i = 0; i < _n; ++ i) {
        double l = *lower_bound(Theta[i].begin(), Theta[i].end(), alpha);
        double r = *(-- upper_bound(Theta[i].begin(), Theta[i].end(), beta));
        if (r < l) result[_data[i].m_id] = 1.0/_obj_cnt[_data[i].m_obj_id];
        else {
            if (rsky_prob[i].count(mp(l, r)) != 0) {
                result[_data[i].m_id] = rsky_prob[i][mp(l, r)];
            }
        }
    }
}


void Dataset::dual_ms_preprocess() {
    _filter_zero_skyprob();
    sort(_data, _data + _n, [](const Instance& t, const Instance& s){
        return t.m_obj_id < s.m_obj_id;
    });
    /* for (int i = 0; i < _n; ++ i) {
        cout << _data[i] << endl;
    } */
    _obj_id_map.resize(_m);
    _id_map.resize(_n);
    int prev = -1, obj_id = -1;
    for (int i = 0; i < _n; ++ i) {
        if (_data[i].m_obj_id != prev) {
            ++ obj_id;
        }
        prev = _data[i].m_obj_id;
        _obj_id_map[obj_id] = prev;
        _data[i].m_obj_id = obj_id;
        _id_map[i] = _data[i].m_id;
        _data[i].m_id = i;
        // cout << _data[i] << endl;
    }
    _m4sp = obj_id + 1;
    // cout << _m4sp << endl;
    _obj_cnt4sp = new int[_m4sp];
    int _real_cnt4sp[_m4sp];
    memset(_real_cnt4sp, 0, _m4sp*sizeof(int));
    for (int i = 0; i < _m4sp; ++ i) {
        _obj_cnt4sp[i] = _obj_cnt[_obj_id_map[i]];
    }
    for (int i = 0; i < _n; ++ i) {
        ++ _real_cnt4sp[_data[i].m_obj_id];
    }
    _kd4sp.resize(_m4sp);
    int start = 0;
    for (int i = 0; i < _m4sp; ++ i) {
        _kd4sp[i] = new KDTree4SP(_dim, _n, _obj_cnt4sp[i]);
        int n = (_n - _real_cnt4sp[i])*_real_cnt4sp[i];
        // cout << "consturct tree for object " << i << " #points: " << n << endl;
        int index = 0;
        ShiftPoint *points = new ShiftPoint[n];
        // cout << n << endl;
        for (int j = start; j < start + _real_cnt4sp[i]; ++ j) {
            // cout << _data[i].m_obj_id << endl;
            for (int k = 0; k < _n; ++ k) {
                if (_data[k].m_obj_id != i) {
                    points[index].m_pt = &_data[j];
                    points[index].m_ps = &_data[k];
                    /* points[index].m_coord.resize(_dim);
                    for (int l = 0; l < _dim; ++ l) {
                        points[index].m_coord[l] = _data[j].m_point[l] - _data[k].m_point[l];
                    } */
                    ++ index;
                }
            }
        }
        start += _real_cnt4sp[i];
        _kd4sp[i]->build(n, points);
    }
    cout << "done\n";
}

// id has been mapped to 1 to _n
void Dataset::dual_ms(HyperRect& ratios, unordered_map<int, double>& result) {
    // unordered_map<int, double> rsky_prob;
    for (int i = 0; i < _n; ++ i) {
        result[i] = 1.0/_obj_cnt4sp[_data[i].m_obj_id];
    }
    for (int i = 0; i < _m4sp; ++ i) {
        _kd4sp[i]->range_query(ratios, result);
    }
    /* for (auto iter : rsky_prob) {
        if (iter.second != 0) {
            result[_id_map[iter.first]] = iter.second;
        }
    } */
}

void Dataset::kdtree_preprocess() {
    _filter_zero_skyprob();
}

void Dataset::kdtree(HyperRect& ratios, unordered_map<int, double>& result) {
    vector<Point> vertices;
    ratios.get_vertices(vertices);
    KDTree kdtree(vertices.size());
    kdtree.load_data(_n, _data, vertices);
    kdtree.mix_sky_prob(_m, _obj_cnt, result);
}