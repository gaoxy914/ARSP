#include "kdtree.h"

void KDTree::_clear(KDTree::Node *node) {
    if (node == nullptr) return;
    if (node->m_left) _clear(node->m_left);
    if (node->m_right) _clear(node->m_right);
    delete node;
}

void KDTree::_mix_sky_prob(int *index, int cnt, int depth, const int*obj_cnt, unordered_map<int, double>& result) {
    if (cnt <= 0) return;
    const int axis = depth%_dim;
    const int mid = cnt/2;
    nth_element(index, index + mid, index + cnt, [&](int lhs, int rhs){
        return _data[lhs].m_point[axis] < _data[rhs].m_point[axis];
    });
    HyperRect mbr(_data[*index].m_point, _data[*index].m_point);
    for (int i = 0; i < cnt; ++ i) {
        mbr.append(_data[index[i]].m_point);
    }
    vector<int> cross_index, dominate_index;
    int obj_id;
    double delta;
    for (int index : _iter.m_index) {
        if (_data[index].dominate(mbr)) {
            dominate_index.push_back(index);
            obj_id = _data[index].m_obj_id;
            delta = obj_cnt[obj_id] - _iter.m_sigma[obj_id];
            ++ _iter.m_sigma[obj_id];
            if (delta > 1) {
                _iter.m_beta *= (delta - 1)/delta;
            } else if (delta == 1) {
                _iter.m_beta *= obj_cnt[obj_id];
                ++ _iter.m_xi;
                break;
            }
        } else if (_data[index].cross(mbr)) {
            cross_index.push_back(index);
        }
    }
    if (_iter.m_xi == 0) {
        _iter.m_index.swap(cross_index);
        if (cnt > 1) {
            _mix_sky_prob(index, mid, depth + 1, obj_cnt, result);
            _mix_sky_prob(index + mid, cnt - mid, depth + 1, obj_cnt, result);
        } else if (cnt == 1) {
            obj_id = _data[*index].m_obj_id;
            double sky_prob = _iter.m_beta/(obj_cnt[obj_id] - _iter.m_sigma[obj_id]);
            result.insert(make_pair(_data[*index].m_id, sky_prob));
        }
        _iter.m_index.swap(cross_index);
    }
    for (int index : dominate_index) {
        obj_id = _data[index].m_obj_id;
        delta = obj_cnt[obj_id] - _iter.m_sigma[obj_id];
        -- _iter.m_sigma[obj_id];
        if (delta > 0) {
            _iter.m_beta *= (delta + 1)/delta;
        } else if (delta == 0) {
            _iter.m_beta /= obj_cnt[obj_id];
            -- _iter.m_xi;
        }
    }
}

void KDTree::_init_iterator(const int& m) {
    _iter.m_sigma = new int[m];
    memset(_iter.m_sigma, 0, m*sizeof(int));
    _iter.m_index.resize(_n);
    iota(begin(_iter.m_index), end(_iter.m_index), 0);
    _iter.m_beta = 1;
    _iter.m_xi = 0;
}

KDTree::KDTree(const int& dim) {
    _dim = dim;
    _n = 0;
    _root = nullptr;
    _data = nullptr;
}

KDTree::~KDTree() {
    if (_root) {
        _clear(_root);
        _root = nullptr;
    }
    if (_iter.m_sigma) {
        delete[] _iter.m_sigma;
        _iter.m_sigma = nullptr;
    }
    if (_data) {
        delete[] _data;
        _data = nullptr;
    }
}

void KDTree::load_data(const int& n, const Instance* data, const vector<Point>& vertices) {
    _n = n;
    _data = new Instance[_n];
    if (vertices.empty()) {
        for (int i = 0; i < _n; ++ i) {
            _data[i] = Instance(data[i]);
        }
    } else {
        for (int i = 0; i < _n; ++ i) {
            _data[i] = Instance(_dim, nullptr, data[i].m_id, data[i].m_obj_id);
            for (int j = 0; j < _dim; ++ j) {
                _data[i].m_point.m_coord[j] = data[i].score(vertices[j]);
            }
        }
    }   
}

void KDTree::mix_sky_prob(const int& m, const int* obj_cnt, unordered_map<int, double>& result) {
    _m = m;
    _init_iterator(m);
    if (_root != nullptr) _clear(_root);
#ifdef _DEBUG_    
    assert(_n > 0 && _data != nullptr);
#endif
    int *index = new int[_n];
    iota(index, index + _n, 0);
    _mix_sky_prob(index, _n, 0, obj_cnt, result);
}