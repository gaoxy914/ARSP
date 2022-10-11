#include "object.hpp"

class DataOperator {
    int dim, m, n, cnt;
    int* ins_cnt;
    vector<InstanceBase*> dataset;
public:
    DataOperator() : dim(0), m(0), n(0), cnt(0) {}

    DataOperator(const int& dim, const int& m) : dim(dim), m(m), n(0), cnt(0) {}

    DataOperator(const int& dim, const int& m, const int& cnt) : dim(dim), m(m), n(0), cnt(cnt) {
        this->ins_cnt = new int[m];
        dataset.resize(m);
    }

    void PrintData() {
        cout << "m = " << m << "\t n = " << n << endl;
        for (int i = 0; i < m; ++ i) {
            cout << "obj " << i << "\t cnt = " << ins_cnt[i] << endl;
            for (int j = 0; j < ins_cnt[i]; ++ j) {
                cout << "id = " << dataset[i][j].ins_id << "\t (";
                for (int k = 0; k < dim - 1; ++ k) cout << dataset[i][j].coord[k] << ", ";
                cout << dataset[i][j].coord[dim - 1] << ")\t";
                cout << "prob = " << 1/dataset[i][j].prob << endl;
            }
        }
    }

    void LoadData(const char *path) {
        ifstream file((string(path) + to_string(dim) + string("/") + to_string(m) + string("/cnt.data")).c_str(), ios::in);
        for (int i = 0; i < m; ++ i) { file >> ins_cnt[i]; n += ins_cnt[i]; }
        file.close();
        file.open((string(path) + to_string(m) + string("/instances.data")).c_str(), ios::in);
        for (int i = 0; i < m; ++ i) {
            dataset[i] = new InstanceBase[ins_cnt[i]];
            for (int j = 0; j < ins_cnt[i]; ++ j) {
                file >> dataset[i][j].obj_id;
                file >> dataset[i][j].ins_id;
                file >> dataset[i][j].prob;
                dataset[i][j].coord = new double[dim];
                for (int k = 0; k < dim; ++ k) file >> dataset[i][j].coord[k];
            }
        }
        file.close();
    }

    void WriteData(const char *path) {
        ofstream file((string(path) + to_string(dim) + string("/") + to_string(m) + string("/cnt.data")).c_str(), ios::out);
        for (int i = 0; i < m; ++ i) file << ins_cnt[i] << " ";
        file.close();
        file.open((string(path) + to_string(dim) + string("/") + to_string(m) + string("/instances.data")).c_str(), ios::out);
        for (int i = 0; i < m; ++ i) {
            for (int j = 0; j < ins_cnt[i]; ++ j) {
                file << dataset[i][j].obj_id << " ";
                file << dataset[i][j].ins_id << " ";
                file << dataset[i][j].prob << " ";
                for (int k = 0; k < dim; ++ k) file << dataset[i][j].coord[k] << " ";
            }
        }
        file.close();
    }

    void GenData(const int& mode, const double& var1 = 0.05, const double& var2 = 0.025) {
        srand((unsigned)time(nullptr));
        switch (mode) {
        case 1:
            for (int i = 0; i < m; ++ i)
                ins_cnt[i] = GenIndePoints(i, dataset[i]);
            break;
        case 2:
            for (int i = 0; i < m; ++ i)
                ins_cnt[i] = GenAntiPoints(i, dataset[i], var1);
            break;
        case 3:
            for (int i = 0; i < m; ++ i)
                ins_cnt[i] = GenCorrPoints(i, dataset[i], var1, var2);
            break;
        default:
            break;
        }
    }
private:
    double RandNormal(const double& med, const double& var) const {
        default_random_engine generator;
        normal_distribution<double> distribution(med, var);
        return distribution(generator);
    }

    /* generate independent instances with uniform distribution */
    int GenIndePoints(const int& obj_id, InstanceBase* &points) {
        double *center = new double[dim];
        for (int i = 0; i < dim; ++ i) center[i] = rand()/double(RAND_MAX);
        int ins_cnt = rand()%cnt + 1;
        double *l = new double[dim];
        double *r = new double[dim];
        for (int i = 0; i < dim; ++ i) {
            double length = RandNormal(1.0, 0.025);
            l[i] = center[i] - length/2 > 0 ? center[i] - length/2 : 0;
            r[i] = center[i] + length/2 < 1 ? center[i] + length/2 : 1;
        }
        points = new InstanceBase[ins_cnt];
        for (int i = 0; i < ins_cnt; ++ i) {
            points[i].obj_id = obj_id;
            // points[i].ins_id = i;
            points[i].ins_id = n ++;
            points[i].coord = new double[dim];
            for (int j = 0; j < dim; ++ j) points[i].coord[j] = rand()/double (RAND_MAX)*(r[j] - l[j]) + l[j];
            // points[i].prob = 1/double(ins_cnt);
            points[i].prob = double(ins_cnt);
        }
        delete[] center;
        delete[] l;
        delete[] r;
        return ins_cnt;
    }

    /* generate anti-correlated instances with uniform distribution */
    int GenAntiPoints(const int& obj_id, InstanceBase* &points, const double& var) {
        double c = 0;
        do {c = RandNormal(0.5, var); } while (c >= 1 || c <= 0);
        double *center = new double[dim];
        do {
            center[dim - 1] = dim*c;
            for (int i = 0; i < dim - 1; ++ i) {
                center[i] = rand()/double(RAND_MAX);
                center[dim - 1] -= center[i];
            }
        } while (center[dim - 1] < 0 || center[dim - 1] > 1);
        int ins_cnt = rand()%cnt + 1;
        double *l = new double[dim];
        double *r = new double[dim];
        for (int i = 0; i < dim; ++ i) {
            double length = RandNormal(1.0, 0.025);
            l[i] = center[i] - length/2 > 0 ? center[i] - length/2 : 0;
            r[i] = center[i] + length/2 < 1 ? center[i] + length/2 : 1;
        }
        points = new InstanceBase[ins_cnt];
        for (int i = 0; i < ins_cnt; ++ i) {
            points[i].obj_id = obj_id;
            // points[i].ins_id = i;
            points[i].ins_id = n ++;
            points[i].coord = new double[dim];
            for (int j = 0; j < dim; ++ j) {
                points[i].coord[j] = rand()/double(RAND_MAX)*(r[j] - l[j]) + l[j];
            }
            // points[i].prob = 1/double(ins_cnt);
            points[i].prob = double(ins_cnt);
        }
        delete[] center;
        delete[] l;
        delete[] r;
        return ins_cnt;
    }

    /* generate correlated instances with unoform distribution */
    int GenCorrPoints(const int& obj_id, InstanceBase* &points, const double& var1, const double& var2) {
        double c = 0;
        do { c = RandNormal(0.5, var1); } while (c >= 1 || c <= 0);
        double *center = new double[dim];
        do {
            center[dim - 1] = dim*c;
            for (int i = 0; i < dim - 1; ++ i) {
                center[i] = RandNormal(0.5, var2);
                center[dim - 1] -= center[i];
            }
        } while (center[dim - 1] < 0 || center[dim - 1] > 1);
        int ins_cnt = rand()%cnt + 1;
        double *l = new double[dim];
        double *r = new double[dim];
        for (int i = 0; i < dim; ++ i) {
            double length = RandNormal(1.0, 0.025);
            l[i] = center[i] - length/2 > 0 ? center[i] - length/2 : 0;
            r[i] = center[i] + length/2 < 1 ? center[i] + length/2 : 1;
        }
        points = new InstanceBase[ins_cnt];
        for (int i = 0; i < ins_cnt; ++ i) {
            points[i].obj_id = obj_id;
            // points[i].ins_id = i;
            points[i].ins_id = n ++;
            points[i].coord = new double[dim];
            for (int j = 0; j < dim; ++ j) {
                points[i].coord[j] = rand()/double(RAND_MAX)*(r[j] - l[j]) + l[j];
            }
            // points[i].prob = 1/double(ins_cnt);
            points[i].prob = double(ins_cnt);
        }
        delete[] center;
        delete[] l;
        delete[] r;
        return ins_cnt;
    }
};

/*
 * ./gendata path data_type dim m cnt var1 var2
 * path: file path
 * data_type:
 * 1 -> independent data
 * 2 -> anti-correlated data (var1 available)
 * 3 -> correlated data (var1, var2 available)
 * dim: dimentionality
 * m: number of uncertain tuples
 * cnt: expected number of instances per uncertain tuples
 */
int main(int argc, char const *argv[]) {
    int mode = atoi(argv[2]);
    int dim = atoi(argv[3]);
    int m = atoi(argv[4]);
    int cnt = atoi(argv[5]);
    DataOperator data_operator(dim, m, cnt);
    double var1 = 0.05, var2 = 0.025;
    if (argc >= 7) var1 = atof(argv[6]);
    if (argc == 8) var2 = atof(argv[7]);
    data_operator.GenData(mode, var1, var2);
    // data_operator.LoadData(argv[1]);
    // data_operator.PrintData();
    data_operator.WriteData(argv[1]);
    return 0;
}
