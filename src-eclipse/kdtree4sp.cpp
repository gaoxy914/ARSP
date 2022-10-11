#include "kdtree4sp.h"

KDTree4SP::Node* KDTree4SP::_build_rec(int *index, int cnt, int depth) {
    if (cnt <= 0) return nullptr;
    const int axis = depth%_dim;
    const int mid = cnt/2;
    nth_element(index, index + mid, index + cnt, [&](int lhs, int rhs){
        return _data[lhs].m_coord[axis] < _data[rhs].m_coord[axis];
    });
    Node *node = new Node();
    node->m_cnt = cnt;
    // node->m_index = new int[cnt];
    node->m_mbr = HyperRect(_dim, _data[*index].m_coord.data(), _data[*index].m_coord.data());
    for (int i = 0; i < cnt; ++ i) {
        // node->m_index[i] = index[i];
        node->m_domiante[_data[index[i]].m_ps->m_id] = true;
        node->m_mbr.append(_data[index[i]].m_coord);
    }
    if (cnt > 1) {
        node->m_left = _build_rec(index, mid, depth + 1);
        node->m_right = _build_rec(index + mid, cnt - mid, depth + 1);
    }
    return node;
}

void KDTree4SP::_range_query_rec(const Node* node, const HyperRect& rect, unordered_map<int, bool>& result) {
    if (node == nullptr) return;
    if (node->m_mbr.upper_edominate(rect)) {
        for (auto iter : node->m_domiante) {
            result[iter.first] = true;
        }
    } else if (node->m_cnt > 1 && node->m_mbr.lower_edominate(rect)) {
        _range_query_rec(node->m_left, rect, result);
        _range_query_rec(node->m_right, rect, result);
    }
}

KDTree4SP::KDTree4SP(const int& dim) {
    _dim = dim;
    _n = 0;
    _root = nullptr;
    _data = nullptr;
    _origin = Instance(_dim);
}

KDTree4SP::~KDTree4SP() {}

void KDTree4SP::build(const int& n, ShiftPoint *data) {
    _n = n;
    _data = data;
    int *index = new int[_n];
    iota(index, index + _n, 0);
    _root = _build_rec(index, _n, 0);
}

void KDTree4SP::range_query(const HyperRect& rect, unordered_map<int, bool>& result) {
    // cout << _root->m_mbr << endl;
    _range_query_rec(_root, rect, result);
}