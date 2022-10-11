#include "eclipse.h"

Dataset::Dataset(const int& dim, const int& n) {
    _dim = dim;
    _n = n;
    _sky = 0;
}

Dataset::~Dataset() {
    if (_data != nullptr) {
        delete[] _data;
    }
}

void Dataset::gen_data() {
    srand((unsigned)time(nullptr));
    srand48((unsigned)time(nullptr));
    _data = new Instance[_n];
    for (int i = 0; i < _n; ++ i) {
        _data[i].m_id = i;
        _data[i].m_point = Point(_dim);
        for (int j = 0; j < _dim; ++ j) {
            _data[i].m_point.m_coord[j] = drand48();
        }
    }
    // compute skyline
    vector<Instance> skyline;
    sort(_data, _data + _n, [&](Instance t, Instance s){
        return t.m_point[0] < s.m_point[0];
    });
    bool in_skyline;
    for (int i = 0; i < _n; ++ i) {
        in_skyline = true;
        for (int j = 0; j < skyline.size(); ++ j) {
            if (skyline[j].dominate(_data[i])) {
                in_skyline = false;
                break;
            }
        }
        if (in_skyline) {
            skyline.push_back(_data[i]);
        }
    }
    // write data
    ofstream file(string("data/inde/") + to_string(_dim) + "_" + to_string(_n) + string(".data").c_str(), ios::out);
    _sky = skyline.size();
    file << _sky << endl;
    for (int i = 0; i < _sky; ++ i) {
        file << _data[i].m_id << " ";
        for (int j = 0; j < _dim - 1; ++ j) {
            file << _data[i].m_point[j] << " ";
        }
        file << _data[i].m_point[_dim - 1] << endl;
    }
}

void Dataset::load_data() {
    ifstream file(string("data/inde/") + to_string(_dim) + "_" + to_string(_n) + string(".data").c_str(), ios::in);
    if (!file.is_open()) {
        printf("Fail in loading data.\n");
        exit(1);
    }
    file >> _sky;
    _data = new Instance[_sky];
    for (int i = 0; i < _sky; ++ i) {
        file >> _data[i].m_id;
        _data[i].m_point = Point(_dim);
        for (int j = 0; j < _dim; ++ j) {
            file >> _data[i].m_point.m_coord[j];
        }
    }
    file.close();
    /* for (int i = 0; i < _sky; ++ i) {
        cout << _data[i] << endl;
    } */
}

void Dataset::build_quadtree() {
    vector<HyperPlane> hyperplanes; // w[1]x[1] + ... + w[d-1]x[d-1] = w[d]
    hyperplanes.resize(_sky*(_sky - 1)/2);
    int index = 0;
    for (int i = 0; i < _sky; ++ i) {
        for (int j = 0; j < i; ++ j) {
            // cout << _data[i].m_id << '\t' << _data[j].m_id << endl;
            hyperplanes[index].m_pt = &_data[i];
            hyperplanes[index].m_ps = &_data[j];
            hyperplanes[index].m_plane.resize(_dim);
            for (int k = 0; k < _dim; ++ k) {
                hyperplanes[index].m_plane[k] = _data[i].m_point[k] - _data[j].m_point[k];
            }
            if (_data[i].m_id == 332 && _data[j].m_id == 236) {
                cout << index << endl;
            }
            ++ index;
        }
    }
    _qd4hp = new QuadTree4HP(_dim - 1);
    // cout << hyperplanes.size() << endl;
    _qd4hp->build(hyperplanes);
}

void Dataset::build_multi_tree() {
    int n = _sky*(_sky - 1);
    ShiftPoint *shiftpoints = new ShiftPoint[n];
    int index = 0;
    for (int i = 0; i < _sky; ++ i) {
        for (int j = 0; j < _sky; ++ j) {
            if (i != j) {
                shiftpoints[index].m_pt = &_data[i];
                shiftpoints[index].m_ps = &_data[j];
                shiftpoints[index].m_coord.resize(_dim);
                for (int k = 0; k < _dim; ++ k) {
                    shiftpoints[index].m_coord[k] = _data[i].m_point[k] - _data[j].m_point[k];
                }
                ++ index;
            }
        }
    }
    _kd4sp = new KDTree4SP(_dim);
    _kd4sp->build(n, shiftpoints);
}

void Dataset::baseline(const HyperRect& rect, vector<int>& result) {
    Point weight(_dim);
    for (int i = 0; i < _dim - 1; ++ i) {
        weight.m_coord[i] = rect.m_lower[i];
    }
    weight.m_coord[_dim - 1] = 1;
    sort(_data, _data + _sky, Instance::Comparator(weight));
    vector<int> eclipse;
    bool in_eclipse;
    for (int i = 0; i < _sky; ++ i) {
        in_eclipse = true;
        for (auto index : eclipse) {
            if (_data[index].edominate(_data[i], rect)) {
                in_eclipse = false;
                break;
            }
        }
        if (in_eclipse) {
            eclipse.push_back(i);
        }
    }
    for (auto index : eclipse) {
        result.push_back(_data[index].m_id);
    }
}

void Dataset::quad(const HyperRect& rect, vector<int>& result) {
    Point weight(_dim);
    weight.m_coord[_dim - 1] = 1;
    for (int i = 0; i < _dim - 1; ++ i) {
        weight.m_coord[i] = rect.m_lower[i];
    }
    sort(_data, _data + _sky, Instance::Comparator(weight));// this sort changes the data of the pointer, so wrong correct not conduct baseline
    unordered_map<int, int> order;
    for (int i = 0; i < _sky; ++ i) {
        order.insert(make_pair(_data[i].m_id, i));
        cout << _data[i].m_id << '\t' << _data[i].score(weight) << '\t' << i << endl;
    }
    result.push_back(_data->m_id);
    // cout << "find order\n";
    vector<HyperPlane> hyperplanes;
    HyperRect range(rect.m_dim);
    for (int i = 0; i < rect.m_dim; ++ i) {
        range.m_lower[i] = -rect.m_upper[i];
        range.m_upper[i] = -rect.m_lower[i];
    }
    // cout << "build query\t" << range << endl;
    _qd4hp->rect_query(range, hyperplanes);
    // cout << hyperplanes.size() << endl;
    for (auto hp : hyperplanes) {
        cout << hp.m_ps->m_id << '\t' << hp.m_pt->m_id << endl;
        if (hp.m_ps->score(weight) < hp.m_pt->score(weight)) {
            -- order[hp.m_pt->m_id];
            if (order[hp.m_pt->m_id] == 0) {
                result.push_back(hp.m_pt->m_id);
            }
        } else {
            -- order[hp.m_ps->m_id];
            if (order[hp.m_ps->m_id] == 0) {
                result.push_back(hp.m_ps->m_id);
            }
        }
    }
}

void Dataset::multi(const HyperRect& rect, unordered_map<int, bool>& result) {
    for (int i = 0; i < _sky; ++ i) {
        result[_data[i].m_id] = false;
    }
    _kd4sp->range_query(rect, result);
}