#ifndef ALL_ECL_PROB_H
#define ALL_ECL_PROB_H

#include "instance.h"
#include "kdtree.h"
#include "kdtree4sp.h"

class Dataset {
    // data
    int _dim; // dimensionality
    int _m; // #uncertain tuple
    int _cnt; // upper bound of #instance of an uncertain tuple
    int _n; // #instance
    double _l; // upper bound of region length
    double _phi_prob;
    string _dpath_prefix;

    int *_obj_cnt;
    Instance *_data;

    vector<KDTree4SP*> _kd4sp;
    vector<int> _obj_id_map;
    vector<int> _id_map;
    int _m4sp;
    int *_obj_cnt4sp;

    vector<vector<double> > Theta;
    vector<map<pair<double, double>, double> > rsky_prob;

    double _roundoff(const double& value, unsigned char prec);
    double _rand_uniform(const double& a, const double& b);
    double _rand_normal(const double& med, const double& var);
    void _gen_inde_data();
    void _gen_anti_data();
    void _gen_corr_data();

    int _filter_zero_skyprob();

public:
    Dataset(const int& dim, const int& m, const int& cnt, const double& l, const double& phi_prob);
    virtual ~Dataset();

    void gen_data(const char* path, const int& type);
    void load_data(const char* path);
    void write_data(const char* path);
    void load_nba_data();

    void dual_ms_preprocess_2d();
    void dual_ms_2d(HyperRect& ratios, unordered_map<int, double>& result);

    void dual_ms_preprocess();
    void dual_ms(HyperRect& ratios, unordered_map<int, double>& result);

    void kdtree_preprocess();
    void kdtree(HyperRect& ratios, unordered_map<int, double>& result);
};


#endif