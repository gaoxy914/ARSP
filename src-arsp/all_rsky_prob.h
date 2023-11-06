#ifndef ALL_RSKY_PROB_H
#define ALL_RSKY_PROB_H

#include "instance.h"
#include "region.h"
#include "kdtree.h"
#include "quadtree.h"
#include "rtree.h"

class Dataset {
    // data
    int _dim; // dimensionality
    int _m; // #uncertain tuple
    int _cnt; // upper bound of #instance of an uncertain tuple
    int _n; // #instance
    double _l; // upper bound of region length
    double _phi_prob;
    string _dpath_prefix;

    int *_start;
    int *_obj_cnt;
    Instance *_data;
    unordered_map<double, double> id_pid;

    RTree *_rtree;

    double _roundoff(const double& value, unsigned char prec);
    double _rand_uniform(const double& a, const double& b);
    double _rand_normal(const double& med, const double& var);
    void _gen_inde_data();
    void _gen_anti_data();
    void _gen_corr_data();
    void _add_empty();
    
    bool _prune(Region& weights, const vector<Instance*>& P, const HyperRect* rect);
    void _update_P(Region& weights, vector<Instance*>& P, Instance* p);

public:
    Dataset(const int& dim, const int& m, const int& cnt, const double& l, const double& phi_prob);
    virtual ~Dataset();

    int get_n();

    void gen_data(const char* path, const int& type);
    void load_data(const char* path);
    void write_data(const char* path);
    void print_data_sketch();
    void load_nba_data();
    void load_car_data();
    void load_iip_data();

    void build_rtree();

    void add_empty(const char *path);
    void gen_data_vary_m(const char* path, const int& type);
    void gen_data_vary_l(const char* path, const int& type);

    void aggregate_rskyline(Region& weights, vector<int>& result);
    void analyze_rsky_prob(Region& weights);
    void kdtree_traverse_star_larger(Region& weights, unordered_map<int, double>& result);

    void enumerate(Region& weights, unordered_map<int, double>& result);
    void enum_rec(Region& weights, vector<int>& pw, int depth, unordered_map<int, int>& cnt);

    void baseline_LP(Region& weights, unordered_map<int, double>& result);
    void baseline_LP_star(Region& weights, unordered_map<int, double>& result);
    void baseline_V(Region& weights, unordered_map<int, double>& result);
    void branch_bound(Region& weights, unordered_map<int, double>& result);
    void branch_bound_trans_on_the_fly(Region& weights, unordered_map<int, double>& result);
    void kdtree_traverse(Region& weights, unordered_map<int, double>& result);
    void kdtree_traverse_star(Region& weights, unordered_map<int, double>& result);
    void quadtree_traverse_star(Region& weights, unordered_map<int, double>& result);
};

#endif