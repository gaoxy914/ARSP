#include "all_rsky_prob.h"

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
        for (int j = start; j < start + _obj_cnt[i]; ++ j) {
            _data[j].m_obj_id = i;
            _data[j].m_id = j;
            _data[j].m_point = Point(_dim);
            for (int k = 0; k < _dim; ++ k) {
                _data[j].m_point.m_coord[k] = _rand_uniform(lower[k], upper[k]);
            }
        }
        start += _obj_cnt[i];
    }
}

void Dataset::_gen_anti_data() {
    double c = 0;
    double center[_dim], lower[_dim], upper[_dim];
    int start = 0;
    double length = 0;
    for (int i = 0; i < _m; ++ i) {
        center[_dim - 1] = 0.5*_dim + _rand_normal(0, 0.05);
        for (int j = 0; j < _dim - 1; ++ j) {
            center[j] = _rand_uniform(0, 1)*min(1.0, center[_dim - 1]);
            center[_dim - 1] -= center[j];
        }
        if (center[_dim - 1] > 1) center[_dim - 1] = 1;
        for (int j = 0; j < _dim; ++ j) {
            do { 
                length = _rand_normal(_l/2.0, _l/8.0);
            } while (length > _l || length <= 0);
            lower[j] = center[j] - length/2 > 0 ? center[j] - length/2 : 0;
            upper[j] = center[j] + length/2 < 1 ? center[j] + length/2 : 1;
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
    } 
}

void Dataset::_gen_corr_data() {
    double c = 0;
    double center[_dim], lower[_dim], upper[_dim];
    int start = 0;
    double length = 0;
    for (int i = 0; i < _m; ++ i) {
        center[0] = _rand_uniform(0, 1);
        for (int j = 1; j < _dim; ++ j) {
            do {
                center[j] = _rand_normal(center[0], 0.05);
            } while (center[j] < 0 || center[j] > 1);
        }
        for (int j = 0; j < _dim; ++ j) {
            do {
                length = _rand_normal(_l/2.0, _l/8.0);
            } while (length > _l || length <= 0);
            lower[j] = center[j] - length/2 > 0 ? center[j] - length/2 : 0;
            upper[j] = center[j] + length/2 < 1 ? center[j] + length/2 : 1;
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
    }
}

void Dataset::_add_empty() {
    vector<int> index(_m);
    iota(index.begin(), index.end(), 0);
    random_shuffle(index.begin(), index.end());
    int pref_sum[_m];
    memset(pref_sum, 0, _m*sizeof(int));
    pref_sum[0] = _obj_cnt[0] - 1;
    for (int i = 1; i < _m; ++ i) {
        pref_sum[i] = pref_sum[i - 1] + _obj_cnt[i];
    }
    for (int i = 0; i < _phi_prob*_m; ++ i) {
        for (int j = 0; j < _dim; ++ j) {
            _data[pref_sum[index[i]]].m_point.m_coord[j] = 1;
        }
    }
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
    _rtree = nullptr;
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

int Dataset::get_n() {
    return _n;
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
    _add_empty();
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

void Dataset::print_data_sketch() {
    printf("Dataset Parameters: dim = %d, m = %d, cnt = %d, l = %.1f,\
     phi = %.1f, n = %d.\n", _dim, _m, _cnt, _l, _phi_prob, _n);
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
        for (int j = 0; j < 9; ++ j) {
            if (j < _dim) file >> _data[i].m_point.m_coord[j];
            else file >> temp;
        }
        id_pid[obj_id] = temp;
    }
    file.close();
    // for (int i = 0; i < _n; ++ i) {
    //     cout << _data[i] << endl;
    // }
}

void Dataset::load_car_data() {
    ifstream file("data/car.cnt", ios::in);
    if (!file.is_open()) {
        printf("Fail in loading data.\n");
        exit(1);
    }
    for (int i = 0; i < _m; ++ i) {
        file >> _obj_cnt[i];
        _n += _obj_cnt[i];
    }
    file.close();
    file.open("data/car.data", ios::in);
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
        for (int j = 0; j < 4; ++ j) {
            if (j < _dim) file >> _data[i].m_point.m_coord[j];
            else file >> temp;
        }
    }
    file.close();
    // for (int i = 0; i < _n; ++ i) {
    //     cout << _data[i] << endl;
    // }
}

void Dataset::load_iip_data() {
    ifstream file("data/iip.cnt", ios::in);
    if (!file.is_open()) {
        printf("Fail in loading data.\n");
        exit(1);
    }
    for (int i = 0; i < _m; ++ i) {
        file >> _obj_cnt[i];
        _n += _obj_cnt[i];
    }
    file.close();
    file.open("data/iip.data", ios::in);
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
        for (int j = 0; j < 2; ++ j) {
            if (j < _dim) file >> _data[i].m_point.m_coord[j];
            else file >> temp;
        }
    }
    file.close();
    // for (int i = 0; i < _n; ++ i) {
    //     cout << _data[i] << endl;
    // }
}


void Dataset::build_rtree() {
    _rtree = new RTree(_dim);
    for (int i = 0; i < _n; ++ i) {
        _rtree->insert(&_data[i]);
    }
}

void Dataset::add_empty(const char *path) {
    vector<int> index(_m);
    iota(index.begin(), index.end(), 0);
    random_shuffle(index.begin(), index.end());
    int pref_sum[_m];
    memset(pref_sum, 0, _m*sizeof(int));
    pref_sum[0] = _obj_cnt[0] - 1;
    for (int i = 1; i < _m; ++ i) {
        pref_sum[i] = pref_sum[i - 1] + _obj_cnt[i];
    }
    for (double p = 0.05; p < 1; p *= 2) {
        for (int i = 0; i < p*_m; ++ i) {
            for (int j = 0; j < _dim; ++ j) {
                _data[pref_sum[index[i]]].m_point.m_coord[j] = 1;
            }
        }
        _dpath_prefix = to_string(_dim) + "_" + to_string(_m) + "_" + to_string(_cnt);
        stringstream stream;
        stream << fixed << setprecision(1) << _l << "_" << fixed << setprecision(2) << p;
        _dpath_prefix += "_" + stream.str();
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
    }
}

void Dataset::gen_data_vary_m(const char* path, const int& type) {
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
    int start = 0;
    vector<int> index;
    for (int m = 2000; m <= 64000; m *= 2) {
        cout << m << '\t' << start << endl;
        vector<int> temp(m - start);
        iota(temp.begin(), temp.end(), start);
        random_shuffle(temp.begin(), temp.end());
        for (int i = 0; i < (m - start)*_phi_prob; ++ i) {
            index.push_back(temp[i]);
        }
        start = m;
        cout << index.size() << endl;
    }
    int pref_sum[_m];
    memset(pref_sum, 0, _m*sizeof(int));
    pref_sum[0] = _obj_cnt[0] - 1;
    for (int i = 1; i < _m; ++ i) {
        pref_sum[i] = pref_sum[i - 1] + _obj_cnt[i];
    }
    for (int i : index) {
        for (int j = 0; j < _dim; ++ j) {
            _data[pref_sum[i]].m_point.m_coord[j] = 1;
        }
    }
    write_data(path);
}

void Dataset::gen_data_vary_l(const char* path, const int& type) {
    for (int i = 0; i < _m; ++ i) {
        _obj_cnt[i] = rand()%_cnt + 1;
        _n += _obj_cnt[i];
    }

    double c = 0;
    double center[_m][_dim];
    double lower[_dim], upper[_dim];
    int start = 0;
    double length = 0;
    
    _data = new Instance[_n];
    switch (type) {
    case 1:
        for (int i = 0; i < _m; ++ i) {
            for (int j = 0; j < _dim; ++ j)
                center[i][j] = _rand_uniform(0, 1);
        }
        break;
    case 2:
        for (int i = 0; i < _m; ++ i) {
            center[i][_dim - 1] = 0.5*_dim + _rand_normal(0.5, 0.05);
            for (int j = 0; j < _dim - 1; ++ j) {
                center[i][j] = _rand_uniform(0, 1)*min(1.0, center[i][_dim - 1]);
                center[i][_dim - 1] -= center[i][j];
            }
            if (center[i][_dim - 1] > 1) center[i][_dim - 1] = 1;
        }
        break;
    case 3:
        for (int i = 0; i < _m; ++ i) {
            center[i][0] = _rand_uniform(0, 1);
            for (int j = 1; j < _dim; ++ j) {
                do {
                    center[i][j] = _rand_normal(center[i][0], 0.05);
                } while (center[i][j] < 0 || center[i][j] > 1);
            }
        }
        break;
    default:
        break;
    }

    for (double l = 0.1; l <= 0.6; l += 0.1) {
        for (int i = 0; i < _m; ++ i) {
            for (int j = 0; j < _dim; ++ j) {
                do { 
                    length = _rand_normal(_l/2.0, _l/8.0);
                } while (length > _l || length <= 0);
                lower[j] = center[i][j] - length/2 > 0 ? center[i][j] - length/2 : 0;
                upper[j] = center[i][j] + length/2 < 1 ? center[i][j] + length/2 : 1;
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
        }
        _dpath_prefix = to_string(_dim) + "_" + to_string(_m) + "_" + to_string(_cnt);
        stringstream stream;
        stream << fixed << setprecision(1) << l << "_" << fixed << setprecision(2) << _phi_prob;
        _dpath_prefix += "_" + stream.str();
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
    }
}

void Dataset::aggregate_rskyline(Region& weights, vector<int>& result) {
    Instance *mean_data = new Instance[_m];
    Instance *var = new Instance[_m];
    int start = 0;
    for (int i = 0; i < _m; ++ i) {
        mean_data[i].m_id = i;
        mean_data[i].m_obj_id = id_pid[i];
        mean_data[i].m_point = Point(_dim);
        for (int j = start; j < start + _obj_cnt[i]; ++ j) {
            for (int k = 0; k < _dim; ++ k) {
                mean_data[i].m_point.m_coord[k] += _data[j].m_point[k];
            }
        }
        var[i].m_point = Point(_dim);
        for (int k = 0; k < _dim; ++ k) {
            mean_data[i].m_point.m_coord[k] /= _obj_cnt[i];
            var[i].m_point.m_coord[k] = 0;
        }
        for (int j = start; j < start + _obj_cnt[i]; ++ j) {
            for (int k = 0; k < _dim; ++ k) {
                var[i].m_point.m_coord[k] += pow(_data[j].m_point[k] - mean_data[i].m_point[k], 2);
            }
        }
        for (int k = 0; k < _dim; ++ k) {
            var[i].m_point.m_coord[k] /= _obj_cnt[i];
        }
        start += _obj_cnt[i];
    }
    cout << "compute aggregate value\n";
    /* for (int i = 0; i < _m; ++ i) {
        cout << mean_data[i] << endl;
        cout << var[i] << endl;
    } */
    vector<int> rsky_index;
    bool in_rsky;
    weights.compute_vertex();
    sort(mean_data, mean_data + _m, Instance::Comparator(weights.get_inner()));
    for (int i = _m - 1; i >= 0; -- i) {
        in_rsky = true;
        for (auto index : rsky_index) {
            if (mean_data[index].fdominate_V_larger(mean_data[i], weights.get_vertex())) {
                in_rsky = false;
                break;
            }
        }
        if (in_rsky == true) {
            rsky_index.push_back(i);
        }
    }
    cout << rsky_index.size() << endl;
    for (auto index : rsky_index) {
        result.push_back(id_pid[mean_data[index].m_id]);
    }
}

void Dataset::analyze_rsky_prob(Region& weights) {

    unordered_map<int, double> result;
    kdtree_traverse_star_larger(weights, result);

    cout << result.size() << endl;

    vector<pair<int, double> > obj_prob;
    vector<vector<pair<Instance, double> > > ins_prob(_m);
    
    int start = 0;
    double sum;
    for (int i = 0; i < _m; ++ i) {
        sum = 0;
        for (int j = start; j < start + _obj_cnt[i]; ++ j) {
            if (result.count(j) != 0) {
                ins_prob[i].push_back(make_pair(_data[j], result[j]));
                sum += result[j];
            }
        }
        start += _obj_cnt[i];
        if (sum != 0) {
            obj_prob.push_back(make_pair(i, sum));
        }
    }

    sort(obj_prob.begin(), obj_prob.end(), [](const pair<int, double>& p1, const pair<int, double>& p2){
        return p1.second > p2.second;
    });

    int obj_id;
    for (int i = 0; i < obj_prob.size(); ++ i) {
        obj_id = obj_prob[i].first;
        cout << obj_id << '\t' << obj_prob[i].second << endl;
        for (int j = 0; j < ins_prob[obj_id].size(); ++ j) {
            cout << ins_prob[obj_id][j].first << '\t' << ins_prob[obj_id][j].second << endl;
        }
    }
}


void Dataset::kdtree_traverse_star_larger(Region& weights, unordered_map<int, double>& result) {
    weights.compute_vertex();
    KDTree kdtree(weights.get_vertex_size());
    kdtree.load_data(_n, _data, &weights);
    unordered_map<int, double> ins_rsky;
    kdtree.mix_sky_prob_larger(_m, _obj_cnt, ins_rsky);
    for (int i = 0; i < _n; ++ i) {
        if (ins_rsky.find(_data[i].m_id) != ins_rsky.end()) {
            result[id_pid[_data[i].m_obj_id]] += ins_rsky[_data[i].m_id];
        }
    }
}

void Dataset::enumerate(Region& weights, unordered_map<int, double>& result) {
    weights.compute_vertex();
    vector<int> pw(_m, 0);
    unordered_map<int, int> cnt;
    _start = new int[_m];
    _start[0] = 0;
    for (int i = 1; i < _m ; ++ i) {
        _start[i] = _start[i - 1] + _obj_cnt[i - 1];
    }
    enum_rec(weights, pw, 0, cnt);
    /* long long nPW = 1;
    for (int i = 0; i < _m; ++ i) {
        nPW *= _obj_cnt[i];
    } */
    for (auto iter : cnt) {
        result[iter.first] = iter.second;
    }
    delete[] _start;
}

void Dataset::enum_rec(Region& weights, vector<int>& pw, int depth, unordered_map<int, int>& cnt) {
    if (depth == _m) {
        Instance *pw_data = new Instance[_m];
        for (int i = 0; i < _m; ++ i) {
            pw_data[i] = _data[pw[i]];
        }
        sort(pw_data, pw_data + _m, Instance::Comparator(weights.get_inner()));
        vector<int> rsky_index;
        bool in_rsky;
        for (int i = 0; i < _m; ++ i) {
            in_rsky = true;
            for (auto index : rsky_index) {
                if (pw_data[index].fdominate_V(pw_data[i], weights.get_vertex())) {
                    in_rsky = false;
                    break;
                }
            }
            if (in_rsky == true) {
                rsky_index.push_back(i);
            }
        }
        for (auto index : rsky_index) {
            ++ cnt[pw_data[index].m_id];
        }
        delete[] pw_data;
    } else if (depth < _m) {
        for (int i = 0; i < _obj_cnt[depth]; ++ i) {
            pw[depth] = _start[depth] + i;
            enum_rec(weights, pw, depth + 1, cnt);
            // pw.pop_back();
        }
    }
}

void Dataset::baseline_LP(Region& weights, unordered_map<int, double>& result) {
    weights.build_lp();
    sort(_data, _data + _n, Instance::Comparator(weights.get_inner()));
    int sigma[_m], obj_id;
    double beta, delta;
    bool full_dominated;
    for (int i = 0; i < _n; ++ i) {
        full_dominated = false;
        memset(sigma, 0, _m*sizeof(int));
        beta = 1;
#ifdef _DEBUG_
        // cout << "id : " << _data[i].m_id << "\t D: \t";
#endif
        for (int j = 0; j < i; ++ j) {
            if (_data[j].m_obj_id != _data[i].m_obj_id && \
            _data[j].fdominate_LP(_data[i], weights.get_lp())) {
#ifdef _DEBUG_
                // cout << _data[j].m_id << '\t';
#endif
                obj_id = _data[j].m_obj_id;
                delta = _obj_cnt[obj_id] - sigma[obj_id];
                if (delta == 1) {
                    full_dominated = true;
                    break;
                } else if (delta > 1) {
                    beta *= (delta - 1)/delta;
                    ++ sigma[obj_id];
                }
            }
        }
#ifdef _DEBUG_
        // cout << endl;
#endif
        if (!full_dominated) {
            result.insert(make_pair(_data[i].m_id, beta/_obj_cnt[_data[i].m_obj_id]));
        }
    }
}

void Dataset::baseline_LP_star(Region& weights, unordered_map<int, double>& result) {
    weights.build_lp();
    sort(_data, _data + _n, Instance::Comparator(weights.get_inner()));
    int sigma[_m], obj_id;
    double beta, delta;
    bool full_dominated;
    Point pivot(weights.get_inner());
    for (int i = 0; i < _n; ++ i) {
        full_dominated = false;
        memset(sigma, 0, _m*sizeof(int));
        beta = 1;
        for (int j = 0; j < i; ++ j) {
            if (_data[j].m_obj_id != _data[i].m_obj_id) {
                if (_data[j].score(pivot) < _data[i].score(pivot)) {
                    if (_data[j].fdominate_LP(_data[i], weights.get_lp())) {
                        obj_id = _data[j].m_obj_id;
                        delta = _obj_cnt[obj_id] - sigma[obj_id];
                        if (delta == 1) {
                            full_dominated = true;
                            break;
                        } else if (delta > 1) {
                            beta *= (delta - 1)/delta;
                            ++ sigma[obj_id];
                        }
                    }
                }
            }
        }
        if (!full_dominated) {
            result.insert(make_pair(_data[i].m_id, beta/_obj_cnt[_data[i].m_obj_id]));
        }
    }   
}

void Dataset::baseline_V(Region& weights, unordered_map<int, double>& result) {
    weights.compute_vertex();
    sort(_data, _data + _n, Instance::Comparator(weights.get_inner()));
    int sigma[_m], obj_id;
    double beta, delta;
    bool full_dominated;
    for (int i = 0; i < _n; ++ i) {
        full_dominated = false;
        memset(sigma, 0, _m*sizeof(int));
        beta = 1;
#ifdef _DEBUG_
        // if (_data[i].m_id == 0) cout << "id : " << _data[i].m_id << "\t D: \t";
#endif
        for (int j = 0; j < i; ++ j) {
            if (_data[j].m_obj_id != _data[i].m_obj_id && \
            _data[j].fdominate_V(_data[i], weights.get_vertex())) {
#ifdef _DEBUG_
                // cout << _data[j] << endl;
#endif
                obj_id = _data[j].m_obj_id;
                delta = _obj_cnt[obj_id] - sigma[obj_id];
                if (delta == 1) {
                    full_dominated = true;
                    break;
                } else if (delta > 1) {
                    beta *= (delta - 1)/delta;
                    ++ sigma[obj_id];
                }
            }
        }
#ifdef _DEBUG_
        // cout << endl;
        if (_data[i].m_id == 0) {
            for (int j = 0; j < _m; ++ j) {
                if (sigma[j] != 0) {
                    cout << "sigma for " << _data[i].m_id << '\t' << j << '\t' << sigma[j] << '\t' << _obj_cnt[j] << endl;
                }
            }
        }
#endif        
        if (!full_dominated) {
            result.insert(make_pair(_data[i].m_id, beta/_obj_cnt[_data[i].m_obj_id]));
        }
    }
}

bool Dataset::_prune(Region& weights, const vector<Instance*>& P, const HyperRect *rect) {
    Instance score_rect = Instance(weights.get_vertex_size());
    for (int i = 0; i < weights.get_vertex_size(); ++ i) {
        score_rect.m_point.m_coord[i] = rect->lower_score(weights.get_vertex(i));
    }
    for (int i = 0; i < P.size(); ++ i) {
        if (P[i]->dominate(score_rect)) {
            return true;
        }
    }
    return false;
}

void Dataset::_update_P(Region& weights, vector<Instance*>& P, Instance* p) {
    bool flag = false;
    for (int i = 0; i < P.size(); ++ i) {
        if (p->fdominate_V(*P[i], weights.get_vertex())) {
            P[i] = P.back();
            P.pop_back();
            -- i;
            flag = true;
        }
    }
    if (flag) P.push_back(p);
}

void Dataset::branch_bound(Region& weights, unordered_map<int, double>& result) {
    weights.compute_vertex();
    typedef pair<RTree::Branch*, int> Node;
    
    struct Comparator {
        Point m_weight;
        Comparator(const Point& weight) : m_weight(weight) {}

        bool operator ()(const Node* t, const Node* s) const {
            return t->first->m_mbr.lower_score(m_weight) > s->first->m_mbr.lower_score(m_weight);
        }
    };

    Comparator comp(weights.get_inner());
    priority_queue<Node *, vector<Node *>, Comparator> heap(comp);
    for (int i = 0; i < _rtree->m_root->m_cnt; ++ i) {
        Node *node = new Node{&_rtree->m_root->m_branch[i], _rtree->m_root->m_level};
        heap.push(node);
    }
    
    vector<RTree*> aggregate_rtree(_m);
    for (int i = 0; i < _m; ++ i) {
        aggregate_rtree[i] = new RTree(_dim);
    }
    vector<Instance*> P;
    P.resize(_m);
    int start = 0;
    for (int i = 0; i < _m; ++ i) {
        P[i] = new Instance(weights.get_vertex_size());
        for (int j = start; j < start + _obj_cnt[i]; ++ j) {
            for (int k = 0; k < weights.get_vertex_size(); ++ k) {
                P[i]->m_point.m_coord[k] = max(P[i]->m_point[k], _data[j].score(weights.get_vertex(k)));
            }
        }
        start += _obj_cnt[i];
    }
    
    while (!heap.empty()) {
        Node *node = heap.top();
        heap.pop();
        if (!_prune(weights, P, &node->first->m_mbr)) {
            if (node->second == 0) {
                Instance *ins = node->first->m_data;
                int obj_id = ins->m_obj_id;
                double prob = 1.0/_obj_cnt[obj_id];
                bool full_dominated = false;
                for (int i = 0; i < _m; ++ i) {
                    if (i != obj_id) {
                        int obj_cnt = aggregate_rtree[i]->range_search(ins, weights.get_vertex());
                        if (obj_cnt == _obj_cnt[i]) {
                            full_dominated = true;
                            break;
                        } else if (obj_cnt < _obj_cnt[i]) {
                            prob *= (double)(_obj_cnt[i] - obj_cnt)/_obj_cnt[i];
                        }
                    }
                }
                if (!full_dominated) {
                    result.insert(make_pair(ins->m_id, prob));
                    aggregate_rtree[obj_id]->aggregate_insert(ins);
                } /* else {
                    _update_P(weights, P, ins);
                } */
            } else if (node->second > 0) {
                for (int i = 0; i < node->first->m_child->m_cnt; ++ i) {
                    if (!_prune(weights, P, &node->first->m_child->m_branch[i].m_mbr)) {
                        Node *next_node = new Node{&node->first->m_child->m_branch[i], node->second - 1};
                        heap.push(next_node);
                    }
                }
            }
        }
    }
}

void Dataset::branch_bound_trans_on_the_fly(Region& weights, unordered_map<int, double>& result) {
    weights.compute_vertex();
    int tdim = weights.get_vertex_size();
    typedef pair<RTree::Branch*, int> Node;
    
    struct Comparator {
        Point m_weight;
        Comparator(const Point& weight) : m_weight(weight) {}

        bool operator ()(const Node* t, const Node* s) const {
            return t->first->m_mbr.lower_score(m_weight) > s->first->m_mbr.lower_score(m_weight);
        }
    };

    Comparator comp(weights.get_inner());
    priority_queue<Node *, vector<Node *>, Comparator> heap(comp);
    for (int i = 0; i < _rtree->m_root->m_cnt; ++ i) {
        Node *node = new Node{&_rtree->m_root->m_branch[i], _rtree->m_root->m_level};
        heap.push(node);
    }

    vector<RTree*> aggregate_rtree(_m);
    for (int i = 0; i < _m; ++ i) {
        aggregate_rtree[i] = new RTree(tdim);
    }

    vector<Instance*> P;
    vector<pair<int, Instance*> > candidate_P;
    candidate_P.reserve(_m);
    for (int i = 0; i < _m; ++ i) {
        Instance *max_corner = new Instance(tdim);
        candidate_P.push_back(make_pair(0, max_corner));
    }

    while (!heap.empty()) {
        Node *node = heap.top();
        heap.pop();
        if (!_prune(weights, P, &node->first->m_mbr)) {
            if (node->second == 0) {
                Instance *ins = node->first->m_data;
                Instance *tins = new Instance(tdim);
                for (int i = 0; i < tdim; ++ i) {
                    tins->m_point.m_coord[i] = ins->score(weights.get_vertex(i));
                }
                int obj_id = ins->m_obj_id;
                double prob = 1.0/_obj_cnt[obj_id];
//                 bool full_dominated = false;
                for (int i = 0; i < _m; ++ i) {
                    if (i != obj_id) {
                        int obj_cnt = aggregate_rtree[i]->range_search(tins);
//                        if (obj_cnt == _obj_cnt[i]) {
//                            full_dominated = true;
//                            break;
//                        } else if (obj_cnt < _obj_cnt[i]) {
                        prob *= (double)(_obj_cnt[i] - obj_cnt)/_obj_cnt[i];
//                        }
                    }
                }
//                if (!full_dominated) {
                result.insert(make_pair(ins->m_id, prob));
                aggregate_rtree[obj_id]->aggregate_insert(tins);
                for (int i = 0; i < tdim; ++ i) {
                    candidate_P[obj_id].second->m_point.m_coord[i] = max(tins->m_point[i], candidate_P[obj_id].second->m_point[i]);
                }
                ++ candidate_P[obj_id].first;
                if (candidate_P[obj_id].first == _obj_cnt[obj_id]) {
                    P.push_back(candidate_P[obj_id].second);
                    // cout << P.size() << endl;
                }
//                }
            } else if (node->second > 0) {
                for (int i = 0; i < node->first->m_child->m_cnt; ++ i) {
                    if (!_prune(weights, P, &node->first->m_child->m_branch[i].m_mbr)) {
                        Node *next_node = new Node{&node->first->m_child->m_branch[i], node->second - 1};
                        heap.push(next_node);
                    }
                }
            }
        }
    }
}

void Dataset::kdtree_traverse(Region& weights, unordered_map<int, double>& result) {
    weights.compute_vertex();
    KDTree kdtree(weights.get_vertex_size());
    kdtree.load_data(_n, _data, &weights);
    kdtree.build();
    kdtree.sky_prob(_m, _obj_cnt, result);
}

void Dataset::kdtree_traverse_star(Region& weights, unordered_map<int, double>& result) {
    weights.compute_vertex();
    KDTree kdtree(weights.get_vertex_size());
    kdtree.load_data(_n, _data, &weights);
    kdtree.mix_sky_prob(_m, _obj_cnt, result);
}

void Dataset::quadtree_traverse_star(Region& weights, unordered_map<int, double>& result) {
    weights.compute_vertex();
    QuadTree quadtree(weights.get_vertex_size());
    quadtree.load_data(_n, _data, &weights);
    quadtree.mix_sky_prob(_m, _obj_cnt, result);
}