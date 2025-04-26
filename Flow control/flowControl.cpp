#include <iostream>
#include <vector>
#include <cassert>
#include <fstream>
#include <iomanip>

double PI = 4 * atan(1);
using namespace std;

class Solution {
    double omega, z0, z00;
    int h, r, q;
    double c12_re = 0, c12_im = -1;
    double y10_re, y10_im;
    double lambda;
    double t, delta;
    vector<pair<double, double>>upper_stage;
    vector<pair<double, double>>lower_stage;
public:
    Solution(int h, int r, int q, double omega, double z0, double z00): 
        h(h), r(r), q(q), omega(omega), z0(z0), z00(z00)
    {
        y10_re = omega * z0;
        y10_im = -z00;
        double square = (y10_re * y10_re + y10_im * y10_im) * (4 * h * h * (r + q) * (r + q) - (y10_re * y10_re + y10_im * y10_im) * omega * omega);
        assert(square >= 0);
        lambda = 2 * h * (r + q) * y10_re * (sqrt(square));
        delta = acos(sqrt(4 * h * h * (r + q) * (r + q) - (y10_re * y10_re + y10_im * y10_im) * omega * omega) / (2 * h * (r + q))) / omega;
        t = acos(y10_im / sqrt((y10_re * y10_re + y10_im * y10_im))) / omega;
        upper_stage.insert(upper_stage.end(), r, make_pair(0, 0));
        lower_stage.insert(lower_stage.end(), q, make_pair(0, 0));
        upper_stage[0] = make_pair(t - delta, t + delta);
        
        for (int l = 1; l < r; l++) {
            double tmp1 = upper_stage[l - 1].first + 2 * PI / omega + 2 * PI * (l + 1) / omega;
            //double tmp2 = upper_stage[l - 1].first + PI / omega + 2 * PI * (l) / omega;
            upper_stage[l] = make_pair(tmp1, tmp1 + 2 * delta);
            //lower_stage[l] = make_pair(tmp2, tmp2 + 2 * delta);
        }
        for (int l = 0; l < q; l++) {
            double t_ = upper_stage[l].first + PI / omega + 2 * PI * (l + 1) / omega;
            lower_stage[l] = make_pair(t_, t_ + 2 * delta);
        }
    }

    double K1() const {
        double k1 = y10_re + h * (r + q) * (c12_im * (cos(omega * upper_stage[0].first) - cos(omega * upper_stage[0].second)) - c12_re * (sin(omega * upper_stage[0].first) - sin(omega * upper_stage[0].second))) / omega;
        return k1;
    }
    double K2() const {
        double k2 = y10_im - h * (r + q) * (c12_re * (cos(omega * upper_stage[0].first) - cos(omega * upper_stage[0].second)) + c12_im * (sin(omega * upper_stage[0].first) - sin(omega * upper_stage[0].second))) / omega;
        return  k2;
    }
    void print_t() const {
        cout << "Upper stages:\n";
        for (auto i : upper_stage) {
            cout << i.first << ' ' << i.second << endl;
        }
        cout << "Lower stages:\n";
        for (auto i : lower_stage) {
            cout << i.first << ' ' << i.second << endl;
        }
    }
    void write_t(string name) const {
        ofstream out;
        out.open(name);
        if (out.is_open()) {
            for (auto i : upper_stage) {
                out << fixed << setprecision(9) << i.first << ' ' << i.second << endl;
            }
            for (auto i : lower_stage) {
                out << i.first << ' ' << i.second << endl;
            }
        }
    }
    void check() const {
        int cnt = 0;
        for (auto i : upper_stage) {
            if (i.first == i.second)
                cnt++;
        }
        cout << "Lower stages:\n";
        for (auto i : lower_stage) {
            if (i.first == i.second)
                cnt++;
        }
        cout << cnt;
    }

    double get_delta() const {
        return delta;
    }

};

int main()
{
    double omega = 10, z0 = 1, z00 = 3, lambda = 0;
    double pi = 4 * atan(1);
    int r = 100, q = 100, h = 3, m = 0;
    Solution res(h, r, q, omega, z0, z00);
    cout << "K1 = " << res.K1() << endl;
    cout << "K2 = " << res.K2() << endl;
    //res.print_t();
    res.write_t("stages.txt");
    cout << res.get_delta();
}
