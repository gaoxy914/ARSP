// #include "object.hpp"
#include "quadtree.hpp"
#include "cuttingtree.hpp"

class Object {
public:
    int obj_id;
    int cnt;
    InstanceBase* instances;

    Object() : obj_id(-1), cnt(0), instances(nullptr) {}
};

class Dataset {
    int m, n;
    int c;
    vector<Object> objects;
    QuadTree quadtree;

public:
    Dataset(const int& m, const int& c) : m(m), n(0), c(c) { objects.resize(m); quadtree = QuadTree(m); }

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

    void DualPreprocessing(const HyperBox& space) {
        // dual to hyperplanes
        vector<HyperPlane> planes(n);
        int j = 0;
        for (auto obj : objects) {
            for (int i = 0; i < obj.cnt; ++ i) {
                planes[j].obj_id = obj.instances[i].obj_id;
                planes[j].ins_id = obj.instances[i].ins_id;
                planes[j].prob = obj.instances[i].prob;
                planes[j].coef = new double[Dim];
                for (int k = 0; k < Dim - 1; ++ k) planes[j].coef[k] = obj.instances[i].coord[k];
                planes[j].coef[Dim - 1] = -obj.instances[i].coord[Dim - 1];
                ++ j;
            }
        }
        // build multi-level tree structure
        quadtree.Build(planes, space);
    }

    void ShiftedDualPreprocessing() {

    }

    void DualAlg(const HyperBox& R, map<int, double>& results) {
        for (auto obj : objects)
            for (int i = 0; i < obj.cnt; ++ i)
                results[obj.instances[i].ins_id] = quadtree.CalProb(R, obj.instances[i]);
        for (auto iter : results)
            if (iter.second != 0)
                cout << "(" << iter.first << ", " << iter.second << ")\n";
    }

    void ShiftedDualAlg(const HyperBox& R, map<int, double>& results) {

    }

};

int main(int argc, char const *argv[]) {
    map<int, double> results;
    Dataset D(10, 4);
    double l[2] = {1, 1};
    double r[2] = {2, 2};
    HyperBox R(2, l, r);
    D.LoadData(argv[1]);
    cout << "finish data load.\n";
    // D.PrintData();
    double s[3] = {0, 0, 0};
    double e[3] = {10, 10, 10};
    HyperBox space(3, s, e);
    D.DualPreprocessing(space);
    cout << "finish preprocessing.\n";
    D.DualAlg(R, results);
#ifdef _LINUX_
    gettimeofday(&end, NULL);
    long long mtime, seconds, useconds;
    seconds = end.tv_sec - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    mtime = seconds*1000000 + useconds;
    printf("Total time is: %lld\n", mtime);
#endif
    return 0;
}