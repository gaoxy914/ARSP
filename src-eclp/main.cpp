#include "all_ecl_prob.h"

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
    if (argc == 8) {
        int dim = atoi(argv[3]);
        int m = atoi(argv[4]);
        int cnt = atoi(argv[5]);
        double l = atof(argv[6]);
        double prob = atof(argv[7]);
        Dataset D(dim, m, cnt, l, prob);
        HyperRect ratio(dim - 1);
        for (int i = 0; i < dim - 1; ++ i) {
            ratio.m_lower[i] = 0.5;
            ratio.m_upper[i] = 1.6;
        }
        unordered_map<int, double> result;
        struct timeval start, end;
        long long mtime, seconds, useconds;
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
        } else if (strcmp(argv[1], "-dual-ms") == 0) {
            D.load_data(argv[2]);
            if (dim > 2) {
                gettimeofday(&start, nullptr);
                D.dual_ms_preprocess();
                gettimeofday(&end, nullptr);
                seconds = end.tv_sec - start.tv_sec;
                useconds = end.tv_usec - start.tv_usec;
                mtime = seconds*1000000 + useconds;
                cout << "preprocessing times: " << mtime << endl;    
                gettimeofday(&start, nullptr);
                D.dual_ms(ratio, result);
                gettimeofday(&end, nullptr);
                seconds = end.tv_sec - start.tv_sec;
                useconds = end.tv_usec - start.tv_usec;
                mtime = seconds*1000000 + useconds;
                cout << "times: " << mtime << endl;
            } else if (dim == 2) {
                gettimeofday(&start, nullptr);
                D.dual_ms_preprocess_2d();
                gettimeofday(&end, nullptr);
                seconds = end.tv_sec - start.tv_sec;
                useconds = end.tv_usec - start.tv_usec;
                mtime = seconds*1000000 + useconds;
                cout << "preprocessing times: " << mtime << endl;    
                gettimeofday(&start, nullptr);
                D.dual_ms_2d(ratio, result);
                gettimeofday(&end, nullptr);
                seconds = end.tv_sec - start.tv_sec;
                useconds = end.tv_usec - start.tv_usec;
                mtime = seconds*1000000 + useconds;
                cout << "times: " << mtime << endl;
            }
            /* for (auto iter : result) {
                cout << iter.first << '\t' << iter.second << endl;
            } */
        } else if (strcmp(argv[1], "-kdtree-star") == 0) {
            D.load_data(argv[2]);
            gettimeofday(&start, nullptr);
            D.kdtree_preprocess();
            gettimeofday(&end, nullptr);
            seconds = end.tv_sec - start.tv_sec;
            useconds = end.tv_usec - start.tv_usec;
            mtime = seconds*1000000 + useconds;
            cout << "preprocessing times: " << mtime << endl;    
            gettimeofday(&start, nullptr);
            D.kdtree(ratio, result);
            gettimeofday(&end, nullptr);
            seconds = end.tv_sec - start.tv_sec;
            useconds = end.tv_usec - start.tv_usec;
            mtime = seconds*1000000 + useconds;
            cout << "times: " << mtime << endl;    
            /* for (auto iter : result) {
                cout << iter.first << '\t' << iter.second << endl;
            } */
        }
    } else {
        int dim = 3;
        HyperRect ratio(dim - 1);
        for (int i = 0; i < dim - 1; ++ i) {
            ratio.m_lower[i] = 0.36;
            ratio.m_upper[i] = 2.75;
        }
        struct timeval start, end;
        long long mtime, seconds, useconds;
        if (strcmp(argv[1], "-nba") == 0) {
            for (int i = 500; i < 3700; i += 500) {
                unordered_map<int, double> result;
                Dataset nba(dim, i, 0, 0, 0);
                nba.load_nba_data();
                gettimeofday(&start, nullptr);
                // nba.kdtree_preprocess();
                // nba.dual_ms_preprocess_2d();
                nba.dual_ms_preprocess();
                gettimeofday(&end, nullptr);
                seconds = end.tv_sec - start.tv_sec;
                useconds = end.tv_usec - start.tv_usec;
                mtime = seconds*1000000 + useconds;
                cout << "m = " << i << "\tpreprocessing times: " << mtime << endl;    
                gettimeofday(&start, nullptr);
                // nba.kdtree(ratio, result);
                // nba.dual_ms_2d(ratio, result);
                nba.dual_ms(ratio, result);
                gettimeofday(&end, nullptr);
                seconds = end.tv_sec - start.tv_sec;
                useconds = end.tv_usec - start.tv_usec;
                mtime = seconds*1000000 + useconds;
                cout << "m = " << i << " times: " << mtime << " result: " << result.size() << endl;
            } 
        }
        /* for (int dim = 2; dim <= 5; ++ dim) {
            HyperRect ratio(dim - 1);
            for (int i = 0; i < dim - 1; ++ i) {
                ratio.m_lower[i] = 0.36;
                ratio.m_upper[i] = 2.75;
            }
            struct timeval start, end;
            long long mtime, seconds, useconds;
            unordered_map<int, double> result;
            Dataset nba(dim, 1000, 0, 0, 0);
            nba.load_nba_data();
            gettimeofday(&start, nullptr);
            // nba.kdtree_preprocess();
            // nba.dual_ms_preprocess_2d();
            nba.dual_ms_preprocess();
            gettimeofday(&end, nullptr);
            seconds = end.tv_sec - start.tv_sec;
            useconds = end.tv_usec - start.tv_usec;
            mtime = seconds*1000000 + useconds;
            cout << "dim = " << dim << "\tpreprocessing times: " << mtime << endl;    
            gettimeofday(&start, nullptr);
            // nba.kdtree(ratio, result);
            // nba.dual_ms_2d(ratio, result);
            nba.dual_ms(ratio, result);
            gettimeofday(&end, nullptr);
            seconds = end.tv_sec - start.tv_sec;
            useconds = end.tv_usec - start.tv_usec;
            mtime = seconds*1000000 + useconds;
            cout << "dim = " << dim << " times: " << mtime << " result: " << result.size() << endl;
        } */
    }
    return 0;
}
