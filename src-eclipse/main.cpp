#include "eclipse.h"

int main(int argc, char const *argv[]) {
    int dim = atoi(argv[1]);
    int n = atoi(argv[2]);
    Dataset D(dim, n);
    // D.gen_data();
    D.load_data();
    // int dim = 3;
    HyperRect ratio(dim - 1);
    for (int i = 0; i < dim - 1; ++ i) {
        ratio.m_lower[i] = 0.36;
        ratio.m_upper[i] = 2.75;
    }
    /* for (int i = 1024; i <= 1048576; i *= 4) {
        Dataset D(dim, i);
        D.load_data();
        struct timeval start, end;
        long long mtime, seconds, useconds;
        
        vector<int> result;
        D.baseline(ratio, result);
        for (auto id : result) {
            cout << id << '\t';
        }
        cout << endl;
        result.clear();
        vector<int>().swap(result);

        D.build_quadtree();
        gettimeofday(&start, nullptr);
        D.quad(ratio, result);
        gettimeofday(&end, nullptr);
        seconds = end.tv_sec - start.tv_sec;
        useconds = end.tv_usec - start.tv_usec;
        mtime = seconds*1000000 + useconds;
        cout << "times: " << mtime << endl;    
        for (auto id : result) {
            cout << id << '\t';
        }
        cout << endl;

        unordered_map<int, bool> dominated;
        D.build_multi_tree();
        // cout << "done built tree\n";
        gettimeofday(&start, nullptr);
        D.multi(ratio, dominated);
        gettimeofday(&end, nullptr);
        seconds = end.tv_sec - start.tv_sec;
        useconds = end.tv_usec - start.tv_usec;
        mtime = seconds*1000000 + useconds;
        cout << "times: " << mtime << endl;    
        for (auto iter : dominated) {
            if (!iter.second) cout << iter.first << '\t';
        }
        cout << endl;
    } */
    
    vector<int> result;
    /* D.baseline(ratio, result);
    for (auto id : result) {
        cout << id << '\t';
    }
    cout << endl; */
    
    struct timeval start, end;
    long long mtime, seconds, useconds;
    
    result.clear();
    D.build_quadtree();
    gettimeofday(&start, nullptr);
    D.quad(ratio, result);
    gettimeofday(&end, nullptr);
    seconds = end.tv_sec - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    mtime = seconds*1000000 + useconds;
    cout << "times: " << mtime << endl;    
    for (auto id : result) {
        cout << id << '\t';
    }
    cout << endl;

    unordered_map<int, bool> dominated;
    D.build_multi_tree();
    // cout << "done built tree\n";
    gettimeofday(&start, nullptr);
    D.multi(ratio, dominated);
    gettimeofday(&end, nullptr);
    seconds = end.tv_sec - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    mtime = seconds*1000000 + useconds;
    cout << "times: " << mtime << endl;    
    for (auto iter : dominated) {
        if (!iter.second) cout << iter.first << '\t';
    }
    cout << endl;
}