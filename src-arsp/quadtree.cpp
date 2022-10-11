#include "quadtree.h"

void QuadTree::_init_iterator(const int& m) {
    _iter.m_sigma = new int[m];
    memset(_iter.m_sigma, 0, m*sizeof(int));
    _iter.m_index.resize(_n);
    iota(begin(_iter.m_index), end(_iter.m_index), 0);
    _iter.m_beta = 1;
    _iter.m_xi = 0;
}

void QuadTree::_clear(Node *node) {
    if (node == nullptr) return;
    /* if (node->m_children != nullptr) {
        for (int i = 0; i < _nchild; ++ i) {
            _clear(node->m_children[i]);
        }
    } */
    delete node;
}

void QuadTree::_mix_sky_prob(vector<int> index, const int& depth, const int* obj_cnt, unordered_map<int, double>& result) {
    if (index.empty()) return;
    // Node *node = new Node();
    // node->m_index = index;
    HyperRect mbr(_data[index.front()].m_point, _data[index.front()].m_point);
    for (int i = 1; i < index.size(); ++ i) {
        mbr.append(_data[index[i]].m_point);
    }
#ifdef _DEBUG_
    cout << depth << '\t' << node->m_mbr << endl;
    for (int i : index) {
        cout << _data[i] << endl;
    }
#endif
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
            }
        } else if (_data[index].cross(mbr)) {
            cross_index.push_back(index);
        }
    }
    _iter.m_index.swap(cross_index);
    if (index.size() == 1 && _iter.m_xi == 0) {
        int id = index.back();
        obj_id = _data[id].m_obj_id;
        double sky_prob = _iter.m_beta/(obj_cnt[obj_id] - _iter.m_sigma[obj_id]);
        result.insert(make_pair(_data[id].m_id, sky_prob));
    } else if (index.size() > 1 && _iter.m_xi == 0) {
        Point mid_point = mbr.get_mid();
        vector<vector<int> > sub_index;
        sub_index.resize(_nchild);
        int k;
        for (int i : index) {
            k = 0;
            for (int j = 0; j < _dim; ++ j) {
                if (_data[i].m_point[j] > mid_point[j]) {
                    k += (1<<j);
                }
            }
            sub_index[k].push_back(i);
        }
        // node->m_children = new Node*[_nchild];
        for (int i = 0; i < _nchild; ++ i) {
            _mix_sky_prob(sub_index[i], depth + 1, obj_cnt, result);
        }
    }
    _iter.m_index.swap(cross_index);
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

/* QuadTree::Node* QuadTree::_build(vector<int> index, const int& depth) {
    if (index.empty()) return nullptr;
    Node *node = new Node();
    node->m_index = index;
    node->m_mbr = HyperRect(_data[index.front()].m_point, _data[index.front()].m_point);
    for (int i = 1; i < index.size(); ++ i) {
        node->m_mbr.append(_data[index[i]].m_point);
    }
#ifdef _DEBUG_
    cout << depth << '\t' << node->m_mbr << endl;
    for (int i : index) {
        cout << _data[i] << endl;
    }
#endif
    if (index.size() > 1) {
        vector<HyperRect> sub_rect;
        Point mid_point(_dim);
        // TODO: only compute mid point and partition
        node->m_mbr.z_split(sub_rect, mid_point);
        vector<vector<int> > sub_index;
        sub_index.resize(sub_rect.size());
        int k;
        for (int i : index) {
            k = 0;
            for (int j = 0; j < _dim; ++ j) {
                if (_data[i].m_point[j] > mid_point[j]) {
                    k += (1<<j);
                }
            }
            sub_index[k].push_back(i);
        }
        node->m_children = new Node*[_nchild];
        for (int i = 0; i < _nchild; ++ i) {
            node->m_children[i] = _build(sub_index[i], depth + 1);
        }
    }
    return node;
}

void QuadTree::_sky_porb(const Node* node, const int* obj_cnt, unordered_map<int, double>& result) {
    vector<int> cross_index, dominate_index;
    int obj_id;
    double delta;
    for (int index : _iter.m_index) {
        if (_data[index].dominate(node->m_mbr)) {
            dominate_index.push_back(index);
            obj_id = _data[index].m_obj_id;
            delta = obj_cnt[obj_id] - _iter.m_sigma[obj_id];
            ++ _iter.m_sigma[obj_id];
            if (delta > 1) {
                _iter.m_beta *= (delta - 1)/delta;
            } else if (delta == 1) {
                _iter.m_beta *= obj_cnt[obj_id];
                ++ _iter.m_xi;
            }
        } else if (_data[index].cross(node->m_mbr)) {
            cross_index.push_back(index);
        }
    }
    _iter.m_index.swap(cross_index);
    if (node->m_index.size() == 1) {
        if (_iter.m_xi == 0) {
            int index = node->m_index.back();
            obj_id = _data[index].m_obj_id;
            double sky_prob = _iter.m_beta/(obj_cnt[obj_id] - _iter.m_sigma[obj_id]);
            result.insert(make_pair(_data[index].m_id, sky_prob));
#ifdef _DEBUG_
            if (_data[index].m_id == 0) {
                for (int i = 0; i < 10; ++ i) {
                    if (_iter.m_sigma[i] != 0) {
                        cout << "sigma for " << _data[index].m_id << '\t' << i << '\t' << _iter.m_sigma[i] << '\t' << obj_cnt[i] << endl;
                    }
                }
                for (int i : _iter.tmp) {
                    cout << _data[i] << endl;
                }
            }
#endif
        }
    } else if (node->m_index.size() > 1) {
        for (int i = 0; i < _nchild; ++ i) {
            if (node->m_children[i] != nullptr && _iter.m_xi == 0) {
                _sky_porb(node->m_children[i], obj_cnt, result);
            }
        } 
    }
    _iter.m_index.swap(cross_index);
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
} */

QuadTree::QuadTree(const int& dim) {
    _dim = dim;
    _n = 0;
    _nchild = pow(2, _dim);
    _root = nullptr;
}

QuadTree::~QuadTree() {
    if (_root) {
        _clear(_root);
        _root = nullptr;
    }
    if (_iter.m_sigma) {
        delete[] _iter.m_sigma;
        _iter.m_sigma = nullptr;
    }
}

void QuadTree::load_data(const int& n, const Instance* data, Region* weights) {
    _n = n;
    _data.resize(_n);
    bool flag;
    int index = 0;
    if (weights) {
        for (int i = 0; i < _n; ++ i) {
            _data[index] = Instance(_dim, nullptr, data[i].m_id, data[i].m_obj_id);
            flag = false;
            for (int j = 0; j < _dim; ++ j) {
                _data[index].m_point.m_coord[j] = data[i].score(weights->get_vertex(j));
                if (_data[index].m_point[j] < 1) {
                    flag = true;
                }
            }
            if (flag) {
                ++ index;
            }
        }
        _n = index;
    } else {
        for (int i = 0; i < _n; ++ i) {
            _data[i] = Instance(data[i]);
        }
    }
}

/* void QuadTree::build() {
    if (_root != nullptr) _clear(_root);
#ifdef _DEBUG_
    assert(_n > 0 && !_data.empty());
#endif
    vector<int> index(_n);
    iota(begin(index), end(index), 0);
    _root = _build(index, 0);
}

void QuadTree::sky_prob(const int& m, const int* obj_cnt, unordered_map<int, double>& result) {
    _init_iterator(m);
#ifdef _DEBUG_    
    assert(_root);
#endif
    _sky_porb(_root, obj_cnt, result);
} */

void QuadTree::mix_sky_prob(const int& m, const int* obj_cnt, unordered_map<int, double>& result) {
    _init_iterator(m);
    if (_root) _clear(_root);
#ifdef _DEBUG_
    assert(_n > 0 && !_data.empty());
#endif
    vector<int> index(_n);
    iota(begin(index), end(index), 0);
    _mix_sky_prob(index, 0, obj_cnt, result);
}