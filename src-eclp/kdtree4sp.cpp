#include "kdtree4sp.h"

KDTree4SP::Node* KDTree4SP::_build_rec(int *index, int cnt, int depth) {
    // cout << cnt << endl;
    if (cnt <= 0) return nullptr;
    const int axis = depth%_dim;
    const int mid = cnt/2;
    nth_element(index, index + mid, index + cnt, [&](int lhs, int rhs){
        // return _data[lhs].m_coord[axis] < _data[rhs].m_coord[axis];
        return _data[lhs][axis] < _data[rhs][axis];
    });
    Node *node = new Node();
    node->m_lower = new double[_dim];
    node->m_upper = new double[_dim];
    for (int i = 0; i < _dim; ++ i) {
        node->m_lower[i] = _data[*index][i];
        node->m_upper[i] = _data[*index][i];
    }
    // node->m_mbr = HyperRect(_dim, _data[*index].m_coord.data(), _data[*index].m_coord.data());
    int s;
    unordered_map<int, int> ins_cnt;
    for (int i = 0; i < cnt; ++ i) {
        s = _data[index[i]].m_ps->m_id;
        if (ins_cnt.count(s) == 0) {
            ins_cnt[s] = 1;
        } else {
            ++ ins_cnt[s];
        }
        for (int j = 0; j < _dim; ++ j) {
            node->m_lower[j] = min(node->m_lower[j], _data[index[i]][j]);
            node->m_upper[j] = max(node->m_upper[j], _data[index[i]][j]);
        }
        // node->m_mbr.append(_data[index[i]].m_coord);
    }
    node->m_sigma.reserve(ins_cnt.size());
    for (auto iter : ins_cnt) {
        node->m_sigma.push_back(make_pair(iter.first, iter.second));
    }
    if (cnt > 1) {
        node->m_left = _build_rec(index, mid, depth + 1);
        node->m_right = _build_rec(index + mid, cnt - mid, depth + 1);
    }
    return node;
}

// point F-dominates origin
bool KDTree4SP::_dominate(double *point, const HyperRect& rect) {
    bool unique = false;
    double min_value = -point[_dim - 1];
    for (int i = 0; i < _dim - 1; ++ i) {
        if (point[i] < 0) {
            unique = true;
            min_value -= point[i]*rect.m_lower[i];
        } else if (point[i] > 0) {
            unique = true;
            min_value -= point[i]*rect.m_upper[i];
        }
    }
    return min_value >= 0 && unique;
}

void KDTree4SP::_range_query_rec(const Node* node, const HyperRect& rect, unordered_map<int, double>& result) {
    if (node == nullptr) return;
    if (_dominate(node->m_upper, rect)) {
    // if (node->m_mbr.upper_edominate(rect)) {
        double delta;
        for (auto iter : node->m_sigma) {
            // if (result[iter.first] != 0) {
            if (result.count(iter.first) != 0) {
                delta = _obj_cnt - _sigma[iter.first];
                if (delta == iter.second) {
                    // result[iter.first] = 0;
                    result.erase(iter.first);
                } else if (delta > iter.second) {
                    result[iter.first] *= (delta - iter.second)/delta;
                    _sigma[iter.first] += iter.second;
                }
            }
        }
    } else if (_dominate(node->m_lower, rect)) {
    // } else if (node->m_mbr.lower_edominate(rect)) {
        _range_query_rec(node->m_left, rect, result);
        _range_query_rec(node->m_right, rect, result);
    }
}

KDTree4SP::KDTree4SP(const int& dim, const int& ori_n, const int& obj_cnt) {
    _dim = dim;
    _n = 0;
    _obj_cnt = obj_cnt;
    _root = nullptr;
    _data = nullptr;
    _origin = Instance(_dim);
    _ori_n = ori_n;
    _sigma = new int[_ori_n];
}

KDTree4SP::~KDTree4SP() {}

void KDTree4SP::build(const int& n, ShiftPoint *data) {
    _n = n;
    _data = data;
    int *index = new int[_n];
    iota(index, index + _n, 0);
    _root = _build_rec(index, _n, 0);
}

void KDTree4SP::range_query(const HyperRect& rect, unordered_map<int, double>& result) {
    // cout << _obj_cnt << endl;
    memset(_sigma, 0, _ori_n*sizeof(int));
    _range_query_rec(_root, rect, result);
}