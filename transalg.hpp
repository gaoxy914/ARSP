#include "object.hpp"
#include "kdtree.hpp"

class Object {
public:
    int obj_id;
    int cnt;
    InstanceBase* instances;
    
    Object() : obj_id(-1), cnt(0), instances(nullptr) {}
};

class Dataset {
    int m, n;
    vector<Object> objects;

public:
    Dataset(const int& m) {
        n = 0;
        this->m = m;
        objects.resize(m);
    }


    void LoadData(const char* path) {
        ifstream file((string(path) + to_string(m) + string("/cnt.data")).c_str(), ios::in);
        for (int i = 0; i < m; ++ i) { file >> objects[i].cnt; n += objects[i].cnt; }
        file.close();
        file.open((string(path) + to_string(m) + string("/instances.data")).c_str(), ios::in);
        for (int i = 0; i < m; ++ i) {
            objects[i].instances = new InstanceBase[objects[i].cnt];
            for (int j = 0; j < objects[i].cnt; ++ j) {
                file >> objects[i].instances[j].obj_id;
                file >> objects[i].instances[j].ins_id;
                file >> objects[i].instances[j].prob;
                objects[i].instances[j].coord = new double[Dim];
                for (int k = 0; k < Dim; ++ k) file >> objects[i].instances[j].coord[k];
            }
        }
        file.close();
    }

    void PrintData() {
        cout << "m = " << m << "\t n = " << n << endl;
        for (int i = 0; i < m; ++ i) {
            cout << "obj " << i << "\t cnt = " << objects[i].cnt << endl;
            for (int j = 0; j < objects[i].cnt; ++ j) {
                cout << "id = " << objects[i].instances[j].ins_id << "\t (";
                for (int k = 0; k < Dim - 1; ++ k) cout << objects[i].instances[j].coord[k] << ", ";
                cout << objects[i].instances[j].coord[Dim - 1] << ")\t";
                cout << "prob = " << 1/objects[i].instances[j].prob << endl;
            }
        }
    }

    void TransAlg(const HyperBox& R, map<int, double>& results) {
        vector<double*> vertices;
        R.GetVertices(vertices);
        int dim = vertices.size();
        vector<InstanceBase> new_instances;
        new_instances.reserve(n);
        for (auto obj : objects) {
            for (int i = 0; i < obj.cnt; ++ i) {
                InstanceBase cur_ins = obj.instances[i], new_ins;
                new_ins.obj_id = cur_ins.obj_id;
                new_ins.ins_id = cur_ins.ins_id;
                new_ins.prob = cur_ins.prob;
                new_ins.coord = new double[dim];
                for (int j = 0; j < dim; ++ j) new_ins.coord[j] = dot(cur_ins.coord, vertices[j]);
                new_instances.push_back(new_ins);
            }
        }
        // build kd-tree
        KDTree kdtree(dim, m, new_instances);
        // cout << "built\n";
        // calculate skyline probability
        results = kdtree.CalSkyPorb();
        // cout << "done\n";
        for (auto iter : results)
            if (iter.second != 0)
                cout << "(" << iter.first << ", " << iter.second << ")\n";
    }
};

/*
 * ./trans path m
 * path: file path
 * m: number of uncertain tuples
 */
int main(int argc, char const *argv[]) {
    map<int, double> results;
    int m = atoi(argv[2]);
    Dataset D(m);
    double l[Dim - 1] = {1, 1};
    double r[Dim - 1] = {2, 2};
    HyperBox R(Dim - 1, l, r);
    D.LoadData(argv[1]);
    // D.PrintData();
#ifdef _LINUX_
    struct timeval start, end;
    gettimeofday(&start, NULL);
#endif
    D.TransAlg(R, results);
#ifdef _LINUX_
    long long mtime, seconds, useconds;
    gettimeofday(&end, NULL);
    seconds = end.tv_sec - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    mtime = seconds*1000000 + useconds;
    printf("Total time is: %lld\n", mtime);
#endif
    return 0;
}
