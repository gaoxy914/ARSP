#include "rtree.h"

void RTree::_get_branches(Node *node, Branch *branch) {
    assert(node->m_cnt == MAXNODES); // only for full nodes
    for (int i = 0; i < MAXNODES; ++ i) {
        _var.m_branch[i] = node->m_branch[i];
    }
    _var.m_branch[MAXNODES] = *branch;
    _var.m_total_mbr = branch->m_mbr;
    for (int i = 0; i < MAXNODES + 1; ++ i) {
        _var.m_total_mbr.append(_var.m_branch[i].m_mbr);
    }
    _var.m_total_vol = _var.m_total_mbr.spherical_volume();
    node->m_level = -1;
    node->m_cnt = 0;
}

bool RTree::_add_branch(Branch *branch, Node *node, Node **new_node) {
    if (node->m_cnt < MAXNODES) {
        node->m_branch[node->m_cnt] = *branch;
        ++ node->m_cnt;
        return false;
    } else {
        _split_node(node, branch, new_node);
        return true;
    }
}

void RTree::_load_nodes(Node *node1, Node *node2) {
    for (int i = 0; i < MAXNODES + 1; ++ i) {
        if (_var.m_partition[i] == 0) {
            _add_branch(&_var.m_branch[i], node1, nullptr);
        } else if (_var.m_partition[i] == 1) {
            _add_branch(&_var.m_branch[i], node2, nullptr);
        } else {
            return;
        }
    }
}

void RTree::_classify(const int& index, const int& group) {
    _var.m_partition[index] = group;
    if (_var.m_cnt[group] == 0) {
        _var.m_mbr[group] = _var.m_branch[index].m_mbr;
    } else {
        _var.m_mbr[group].append(_var.m_branch[index].m_mbr);
    }
    _var.m_vol[group] = _var.m_mbr[group].spherical_volume();
    ++ _var.m_cnt[group];
}

void RTree::_pick_seeds() {
    int seed0, seed1;
    double worst, waste;
    double volume[MAXNODES + 1];
    HyperRect rect;
    for (int i = 0; i < MAXNODES + 1; ++ i) {
        volume[i] = _var.m_branch[i].m_mbr.spherical_volume();
    }
    worst = -1*_var.m_total_vol - 1;
    for (int i = 0; i < MAXNODES; ++ i) {
        for (int j = i + 1; j < MAXNODES + 1; ++ j) {
            rect = _var.m_branch[i].m_mbr + _var.m_branch[j].m_mbr;
            waste = rect.spherical_volume() - volume[i] - volume[j];
            if (waste > worst) {
                worst = waste;
                seed0 = i;
                seed1 = j;
            }
        }
    }
    _classify(seed0, 0);
    _classify(seed1, 1);
}

void RTree::_choose_partition() {
    double biggest_diff, growth0, growth1, diff;
    int group, chosen, better_group;
    HyperRect rect0, rect1;
    for (int i = 0; i < MAXNODES + 1; ++ i) {
        _var.m_partition[i] = -1;
    }
    _var.m_cnt[0] = 0;
    _var.m_cnt[1] = 0;
    _var.m_vol[0] = 0;
    _var.m_vol[1] = 0;
    _pick_seeds();
    while (_var.m_cnt[0] + _var.m_cnt[1] < MAXNODES + 1 && _var.m_cnt[0] < MAXNODES + 1 - MINNODES \
        && _var.m_cnt[1] < MAXNODES + 1 - MINNODES) {
        biggest_diff = -1;
        for (int i = 0; i < MAXNODES + 1; ++ i) {
            if (_var.m_partition[i] == -1) {
                rect0 = _var.m_branch[i].m_mbr + _var.m_mbr[0];
                rect1 = _var.m_branch[i].m_mbr + _var.m_mbr[1];
                growth0 = rect0.spherical_volume() - _var.m_vol[0];
                growth1 = rect1.spherical_volume() - _var.m_vol[1];
                diff = growth1 - growth0;
                if (diff >= 0) { group = 0; }
                else { group = 1; diff = -diff; }
                if (diff > biggest_diff) {
                    biggest_diff = diff;
                    chosen = i;
                    better_group = group;
                } else if (diff == biggest_diff && _var.m_vol[group] < _var.m_vol[better_group]) {
                    chosen = i;
                    better_group = group;
                }
            }
        }
        _classify(chosen, better_group);
    }
    if (_var.m_cnt[0] + _var.m_cnt[1] < MAXNODES + 1) {
        if (_var.m_cnt[0] >= MAXNODES + 1 - MINNODES) {
            group = 1;
        } else {
            group = 0;
        }
        for (int i = 0; i < MAXNODES + 1; ++ i) {
            if (_var.m_partition[i] == -1) {
                _classify(i, group);
            }
        }
    }
}

void RTree::_split_node(Node *node, Branch *branch, Node **new_node) {
    int level;
    level = node->m_level;
    _get_branches(node, branch);
    _choose_partition();
    *new_node = new RTree::Node();
    (*new_node)->m_cnt = 0;
    (*new_node)->m_level = node->m_level = level;
    _load_nodes(node, *new_node);

}

HyperRect RTree::_cal_mbr(Node *node) {
    bool first_time = true;
    HyperRect mbr;
    for (int i = 0; i < node->m_cnt; ++ i) {
        if (first_time) {
            mbr = node->m_branch[i].m_mbr;
            first_time = false;
        } else {
            mbr = mbr + node->m_branch[i].m_mbr;
        }
    }
    return mbr;
}

int RTree::_pick_branch(HyperRect *rect, Node *node) {
    bool first_time = true;
    double increase;
    double best_increase = -1;
    double volume, best_volume;
    int best;
    HyperRect temp_rect;
    for (int i = 0; i < node->m_cnt; ++ i) {
        volume = node->m_branch[i].m_mbr.spherical_volume();
        temp_rect = *rect + node->m_branch[i].m_mbr;
        increase = temp_rect.spherical_volume() - volume;
        if (increase < best_increase || first_time) {
            best = i;
            best_volume = volume;
            best_increase = increase;
            first_time = false;
        } else if (increase == best_increase && volume < best_volume) {
            best = i;
            best_volume = volume;
            best_increase = increase;
        }
    }
    return best;
}

bool RTree::_insert_rec(HyperRect *rect, Instance *data, Node *node, Node **new_node, const int& level) {
    int index;
    Branch *branch;
    Node *other_node;
    if (node->m_level > level) {
        index = _pick_branch(rect, node);
        if (!_insert_rec(rect, data, node->m_branch[index].m_child, &other_node, level)) {
            if (_aggregated) ++ node->m_branch[index].m_value;
            node->m_branch[index].m_mbr.append(*rect);
            return false;
        } else {
            node->m_branch[index].m_mbr = _cal_mbr(node->m_branch[index].m_child);
            branch = new Branch();
            branch->m_child = other_node;
            branch->m_mbr = _cal_mbr(other_node);
            if (_aggregated) {
                branch->m_value = 0;
                for (int i = 0; i < other_node->m_cnt; ++ i) {
                    branch->m_value += other_node->m_branch[i].m_value;
                }
                node->m_branch[index].m_value -= branch->m_value;
                ++ node->m_branch[index].m_value;
            }
            return _add_branch(branch, node, new_node);
        }
    } else if (node->m_level == level) {
        branch = new Branch();
        branch->m_mbr = *rect;
        branch->m_data = data;
        if (_aggregated) branch->m_value = 1;
        return _add_branch(branch, node, new_node);
    } else {
        return false;
    }
}

bool RTree::_insert(HyperRect *rect, Instance *data, Node **root, const int& level) {
    Node *new_root;
    Node *new_node;
    Branch *branch;
    if (_insert_rec(rect, data, *root, &new_node, level)) {
        new_root = new Node();
        new_root->m_cnt = 0;
        new_root->m_level = (*root)->m_level + 1;
        branch = new Branch();
        branch->m_mbr = _cal_mbr(*root);
        if (_aggregated) {
            branch->m_value = 0;
            for (int i = 0; i < (*root)->m_cnt; ++ i) {
                branch->m_value += (*root)->m_branch[i].m_value;
            }
        }
        branch->m_child = *root;
        _add_branch(branch, new_root, nullptr);
        branch = new Branch();
        branch->m_mbr = _cal_mbr(new_node);
        if (_aggregated) {
            branch->m_value = 0;
            for (int i = 0; i < new_node->m_cnt; ++ i) {
                branch->m_value += new_node->m_branch[i].m_value;
            }
        }
        branch->m_child = new_node;
        _add_branch(branch, new_root, nullptr);
        *root = new_root;
        return true;
    }
    return false;
}

int RTree::_range_search_rec(Node *node, Instance *data) const {
    int cnt = 0;
    for (int i = 0; i < node->m_cnt; ++ i) {
        if (node->m_branch[i].m_mbr.upper_dominate(data->m_point)) {
            cnt += node->m_branch[i].m_value;
        } else if (node->m_level > 0 && node->m_branch[i].m_mbr.lower_dominate(data->m_point)) {
            cnt += _range_search_rec(node->m_branch[i].m_child, data);
        }
    }
    return cnt;
}

int RTree::_range_search_rec(Node *node, Instance *data, const vector<Point>& weights) const {
    int cnt = 0;
    for (int i = 0; i < node->m_cnt; ++ i) {
        if (node->m_branch[i].m_mbr.upper_fdominate_V(data->m_point, weights)) {
            cnt += node->m_branch[i].m_value;
        } else if (node->m_level > 0 && node->m_branch[i].m_mbr.lower_fdominate_V(data->m_point, weights)) {
            cnt += _range_search_rec(node->m_branch[i].m_child, data, weights);
        }
    }
    return cnt;
}

RTree::RTree(const int& dim) {
    m_dim = dim;
    m_root = new Node();
    m_root->m_level = 0;
    _aggregated = false;
}

RTree::~RTree() {}

void RTree::insert(Instance *data) {
    HyperRect rect(data->m_point, data->m_point);
    _insert(&rect, data, &m_root, 0);
}

void RTree::aggregate_insert(Instance *data) {
    _aggregated = true;
    HyperRect rect(data->m_point, data->m_point);
    _insert(&rect, data, &m_root, 0);
    _aggregated = false;
}

int RTree::range_search(Instance *data) {
    if (m_root->m_cnt == 0) return 0;
    return _range_search_rec(m_root, data);
}

int RTree::range_search(Instance *data, const vector<Point>& weights) const {
    if (m_root->m_cnt == 0) return 0;
    return _range_search_rec(m_root, data, weights);
}

void RTree::show() {
    deque<Node *> q;
    q.push_back(m_root);
    int cur_level = m_root->m_level;
    cout << cur_level << ":\n";
    while (!q.empty()) {
        Node *node = q.front();
        q.pop_front();
        if (node->m_level == 0) {
            if (node->m_level < cur_level) {
                cur_level = node->m_level;
                cout << endl << cur_level << ":\n";
            }
            for (int i = 0; i < node->m_cnt; ++ i) {
                cout << "(" << i << ", " << node->m_branch[i].m_mbr << ")\t";
            }
            cout << endl;
        } else {
            if (node->m_level < cur_level) {
                cur_level = node->m_level;
                cout << endl << cur_level << ":\n";
            }
            for (int i = 0; i < node->m_cnt; ++ i) {
                cout << "(" << node->m_branch[i].m_mbr << ", " << node->m_branch[i].m_value << ")\t";
                q.push_back(node->m_branch[i].m_child);
            }
            cout << endl;
        }
    }
}