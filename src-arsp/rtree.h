#ifndef RTREE_H
#define RTREE_H

#define MAXNODES 8
#define MINNODES MAXNODES / 2

#include "hyperrect.h"
#include "instance.h"

class RTree {
public:
    struct Node;

    struct Branch {
        HyperRect m_mbr;
        // aggregation value
        int m_value;
        union {
            Node *m_child;
            Instance *m_data;
        };
    };
    
    struct BranchComparator {
        Point m_weight;
        BranchComparator(const Point& weight) : m_weight(weight) {}

        bool operator ()(const Branch& t, const Branch& s) const {
            return t.m_mbr.lower_score(m_weight) < s.m_mbr.lower_score(m_weight);
        }
    };

    struct Node {
        int m_level;
        int m_cnt;
        Branch m_branch[MAXNODES];
    };

    int m_dim;
    Node *m_root;
    
protected:
    struct PartitionVars {
        int m_partition[MAXNODES + 1];
        Branch m_branch[MAXNODES + 1];
        HyperRect m_total_mbr;
        double m_total_vol;
        int m_cnt[2];
        HyperRect m_mbr[2];
        double m_vol[2];
    };

    PartitionVars _var;
    bool _aggregated;

    void _get_branches(Node *node, Branch *branch);
    bool _add_branch(Branch *branch, Node *node, Node **new_node);
    void _load_nodes(Node *node1, Node *node2);
    void _classify(const int& index, const int& group);
    void _pick_seeds();
    void _choose_partition();
    void _split_node(Node *node, Branch *branch, Node **new_node);
    HyperRect _cal_mbr(Node *node);
    int _pick_branch(HyperRect *rect, Node *node);
    bool _insert_rec(HyperRect *rect, Instance *data, Node *node, Node **new_node, const int& level);
    bool _insert(HyperRect *rect, Instance *data, Node **root, const int& level);
    int _range_search_rec(Node *node, Instance *data) const;
    int _range_search_rec(Node *node, Instance *data, const vector<Point>& weights) const;
public:
    RTree(const int& dim);
    virtual ~RTree();
    void insert(Instance *data);
    void aggregate_insert(Instance *data);
    int range_search(Instance *data);
    int range_search(Instance *data, const vector<Point>& weights) const;

    void show();
};

#endif