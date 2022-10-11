#include "quadtree4hp.h"

QuadTree4HP::Node* QuadTree4HP::_build_rec(HyperRect& space, vector<int> indices, int depth) {
    if (indices.empty()) return nullptr;
    // cout << "depth: " << depth << endl;
    Node *node = new Node();
    node->_indices.swap(indices);
    if (node->_indices.size() > MAXPLANES) {
        node->_children.resize(_nchild);
        for (int i = 0; i < _nchild; ++ i) {
            HyperRect subspace = space.get_subrect(i);
            vector<int> subindices;
            for (int id : node->_indices) {
                if (subspace.intersect(_hyperplanes[id].m_plane)) {
                    subindices.push_back(id);
                }
            }
            node->_children[i] = _build_rec(subspace, subindices, depth + 1);
        }
    }
    return node;
}

void QuadTree4HP::_query_rec(Node *node, const HyperRect& rect, HyperRect& space, vector<int>& result) {
    if (node == nullptr) return;
    if (rect.contain(space)) {
        // cout << "contain " << space << endl;
        for (int id : node->_indices) {
            if (id == 302) {
                cout << _hyperplanes[id].m_ps->m_id << '\t' << _hyperplanes[id].m_pt->m_id << endl;
                cout << "here\n";
            }
            result.push_back(id);
        }
    } else if (rect.intersect(space)) {
        // cout << "intersects " << space << endl;
        if (node->_children.empty()) { // reach leaf node
            for (int id : node->_indices) {
                if (rect.intersect(_hyperplanes[id].m_plane)) {
                    if (id == 302) {
                        cout << _hyperplanes[id].m_ps->m_id << '\t' << _hyperplanes[id].m_pt->m_id << endl;
                        cout << "here\n";
                    }
                    result.push_back(id);
                }
            }
        } else {
            for (int i = 0; i < _nchild; ++ i) {
                HyperRect subspace = space.get_subrect(i);
                // cout << "child " << i << '\t' << subspace << endl;
                _query_rec(node->_children[i], rect, subspace, result);
            }
        }
    }
}

QuadTree4HP::QuadTree4HP(const int& dim) {
    _dim = dim;
    _nchild = pow(2, _dim);
    _space = HyperRect(_dim);
    for (int i = 0; i < _dim; ++ i) {
        _space.m_lower[i] = -10;
    }
}

QuadTree4HP::~QuadTree4HP() {}

void QuadTree4HP::build(const vector<HyperPlane>& hyperplanes) {
    _hyperplanes = hyperplanes;
    _n = _hyperplanes.size();
    for (int i = 0; i < _n; ++ i) {
        cout << i << '\t' << _hyperplanes[i].m_ps->m_id << '\t' << _hyperplanes[i].m_pt->m_id << endl;
    }
    vector<int> indices(_n);
    iota(indices.begin(), indices.end(), 0);
    _root = _build_rec(_space, indices, 0);
    for (int i = 0; i < _n; ++ i) {
        cout << i << '\t' << _hyperplanes[i].m_ps->m_id << '\t' << _hyperplanes[i].m_pt->m_id << endl;
    }
}

void QuadTree4HP::rect_query(const HyperRect& rect, vector<HyperPlane>& result) {
    cout << "query\n";
    vector<int> indices;
    _query_rec(_root, rect, _space, indices);
    // cout << indices.size() << endl;
    sort(indices.begin(), indices.end());
    indices.erase(unique(indices.begin(), indices.end()), indices.end());
    for (int index : indices) {
        result.push_back(_hyperplanes[index]);
    }
}