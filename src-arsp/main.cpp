#include "all_rsky_prob.h"

#define _EFFICIENCY_

/* run ./main op path dim m cnt l p c
 * op : operator
 * path : data path
 * dim : dimensionality of data
 * m : number of objects
 * cnt : maximum number of instance of each object
 * l : length of possible region of each object
 * p : prob of the existence of an empty instance
 * c : number of constraints
 */
int main(int argc, char const *argv[]) {
    if (argc == 9) {
        int dim = atoi(argv[3]);
        int m = atoi(argv[4]);
        int cnt = atoi(argv[5]);
        double l = atof(argv[6]);
        double prob = atof(argv[7]);
        int c = atoi(argv[8]);
        Dataset D(dim, m, cnt, l, prob);
        Region query(dim, c);
        unordered_map<int, double> result;
        struct timeval start, end;
        long long mtime, seconds, useconds;
        D.load_data(argv[2]);
        // query.load_query();
        query.load_weak_rankings();
        D.build_rtree();
#ifdef _EFFICIENCY_
        gettimeofday(&start, nullptr);
#endif
        if (strcmp(argv[1], "-baseline-lp") == 0) {
            D.baseline_LP(query, result);
        } else if (strcmp(argv[1], "-baseline-lp-star") == 0) {
            D.baseline_LP_star(query, result);
        } else if (strcmp(argv[1], "-baseline-V") == 0) {
            D.baseline_V(query, result);
        } else if (strcmp(argv[1], "-branchbound") == 0) {
            D.branch_bound_trans_on_the_fly(query, result);
        } else if (strcmp(argv[1], "-kdtree") == 0) {
            D.kdtree_traverse(query, result);
        } else if (strcmp(argv[1], "-kdtree-star") == 0) {
            D.kdtree_traverse_star(query, result);
        } else if (strcmp(argv[1], "-quadtree-star") == 0) {
            D.quadtree_traverse_star(query, result);
        } else {
            printf("Option Error.\n");
        }
#ifdef _EFFICIENCY_
        gettimeofday(&end, nullptr);
        seconds = end.tv_sec - start.tv_sec;
        useconds = end.tv_usec - start.tv_usec;
        mtime = seconds*1000 + useconds/1000.0;
        cout << argv[1] << " query time: " << mtime << "ms, size: " << result.size() << endl;
        string record_path = "record/";
        if (strcmp(argv[2], "data/inde/") == 0) {
            record_path += "time-inde.csv";
        } else if (strcmp(argv[2], "data/anti/") == 0) {
            record_path += "time-anti.csv";
        } else if (strcmp(argv[2], "data/corr/") == 0) {
            record_path += "time-corr.csv";
        }
        ofstream file(record_path.c_str(), ios::app);
        file << argv[1] << "," << dim << "," << m << "," << cnt << "," << D.get_n() << "," \
         << l << "," << prob << "," << c << "," << mtime << "," << result.size() << "," << query.get_vertex_size() << endl;
        file.close();
#endif
#ifdef _DEBUG_
        cout << "size: " << result.size() << endl;
        for (auto iter : result) {
            cout << iter.first << '\t' << iter.second << endl;
        }
#endif
    } else if (argc == 8) {
        int dim = atoi(argv[3]);
        int m = atoi(argv[4]);
        int cnt = atoi(argv[5]);
        double l = atof(argv[6]);
        double prob = atof(argv[7]);
        Dataset D(dim, m, cnt, l, prob);
        if (strcmp(argv[1], "-gendata") == 0) {
            int type = 0;
            if (strcmp(argv[2], "data/inde/") == 0) {
                type = 1;
            } else if (strcmp(argv[2], "data/anti/") == 0) {
                type = 2;
            } else if (strcmp(argv[2], "data/corr/") == 0) {
                type = 3;
            }
            D.gen_data(argv[2], type);
        } else if (strcmp(argv[1], "-loaddata") == 0) {
            D.load_data(argv[2]);
        } else if (strcmp(argv[1], "-addempty") == 0) {
            prob = 0;
            D.load_data(argv[2]);
            D.print_data_sketch();
            D.add_empty(argv[2]);
        }
    } else if (argc == 4) {
        int dim = atoi(argv[2]);
        int c = atoi(argv[3]);
        Region query(dim, c);
        // Dataset D(dim, c);
        if (strcmp(argv[1], "-genquery") == 0) {
            // query.gen_query();
            query.gen_weak_rankings();
            query.compute_vertex();
            query.print();
        } else if (strcmp(argv[1], "-loadquery") == 0) {
            // query.load_query();
            query.load_weak_rankings();
            query.compute_vertex();
            query.print();
        } else if (strcmp(argv[1], "-addplane") == 0) {
            query.load_query();
            query.print();
            query.add_plane();
            query.print();
        } else if (strcmp(argv[1], "-gendata") == 0) {
            // D.gen_data();
        } else if (strcmp(argv[1], "-loaddata") == 0) {
            // D.load_data();            
        }
    } else if (argc == 3) {
        if (strcmp(argv[1], "-varym") == 0) {
            Dataset D(4, 64000, 400, 0.2, 0.1);
            int type = 0;
            if (strcmp(argv[2], "data/inde/") == 0) {
                type = 1;
            } else if (strcmp(argv[2], "data/anti/") == 0) {
                type = 2;
            } else if (strcmp(argv[2], "data/corr/") == 0) {
                type = 3;
            }
            D.gen_data_vary_m(argv[2], type);
        }
    } else if (argc == 2) {
        if (strcmp(argv[1], "-nba") == 0) {
            int dim = 5;
            int m = 3707;
            Region query(dim, 4);
            query.load_nba_query();
            // query.compute_vertex();

            Dataset D(dim, m, 0, 0, 0);
            unordered_map<int, double> result;
            D.load_nba_data();
            D.analyze_rsky_prob(query);
            /* vector<int> aggregate_player;
            D.aggregate_rskyline(query, aggregate_player);
            ofstream file("aggregate_rsky.txt", ios::out);
            for (auto i : aggregate_player) {
                file << i << endl;
            }
            file.close(); */
        }
    }
    return 0;
}