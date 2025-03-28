///////// Descrete set for possible Assets and possible Incomes

////////////////////////////////////////////////////// Block 1 //////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace std;

// Utility function
double utility(double c, double sigma) {
    if (c <= 0) return -1e10; // penalize infeasible consumption
    if (sigma == 1.0) return log(c);
    return pow(c, 1.0 - sigma) / (1.0 - sigma);
}


////////////////////////////////////////////////////// Block 2 //////////////////////////////////////////////////////
const int na = 500;               // number of asset grid points
const int ny = 2;                 // number of income states
const double beta = 0.96;         // discount factor
const double sigma = 2.0;         // risk aversion
const double r = 0.04;            // interest rate
const double a_min = 0.0;         // borrowing constraint
const double a_max = 50.0;        // max asset
const double tol = 1e-6;          // tolerance
const int max_iter = 1000;        // max iterations

vector<double> agrid(na);         // asset grid

void create_asset_grid() {
    for (int i = 0; i < na; ++i) {
        agrid[i] = a_min + i * (a_max - a_min) / (na - 1);
    }
}


////////////////////////////////////////////////////// Block 3 //////////////////////////////////////////////////////
vector<double> ygrid = {0.5, 1.5};
vector<vector<double> > P = {{0.9, 0.1}, {0.1, 0.9}};





////////////////////////////////////////////////////// Block 4 //////////////////////////////////////////////////////
vector<vector<double> > V(ny, vector<double>(na, 0.0));
vector<vector<int> > policy_idx(ny, vector<int>(na, 0));

void solve_value_function() {
    for (int it = 0; it < max_iter; ++it) {
        double max_diff = 0.0;
        vector<vector<double> > V_new = V;

        for (int iy = 0; iy < ny; ++iy) {
            for (int ia = 0; ia < na; ++ia) {
                double max_val = -1e10;
                int best_a_idx = 0;

                for (int ja = 0; ja < na; ++ja) {
                    double c = ygrid[iy] + (1 + r) * agrid[ia] - agrid[ja];
                    double u = utility(c, sigma);

                    double EV = 0.0;
                    for (int jy = 0; jy < ny; ++jy) {
                        EV += P[iy][jy] * V[jy][ja];
                    }

                    double val = u + beta * EV;

                    if (val > max_val) {
                        max_val = val;
                        best_a_idx = ja;
                    }
                }

                V_new[iy][ia] = max_val;
                policy_idx[iy][ia] = best_a_idx;

                max_diff = max(max_diff, fabs(V_new[iy][ia] - V[iy][ia]));
            }
        }

        V = V_new;
        if (max_diff < tol) break;
    }
}



////////////////////////////////////////////////////// Block 5 //////////////////////////////////////////////////////
vector<vector<double> > stationary_dist(ny, vector<double>(na, 0.0));

void compute_stationary_distribution() {
    // Initial guess: equal probability
    for (int iy = 0; iy < ny; ++iy)
        for (int ia = 0; ia < na; ++ia)
            stationary_dist[iy][ia] = 1.0 / (ny * na);

    for (int iter = 0; iter < 1000; ++iter) {
        vector<vector<double> > new_dist(ny, vector<double>(na, 0.0));

        for (int iy = 0; iy < ny; ++iy) {
            for (int ia = 0; ia < na; ++ia) {
                int ja = policy_idx[iy][ia];
                for (int jy = 0; jy < ny; ++jy) {
                    new_dist[jy][ja] += stationary_dist[iy][ia] * P[iy][jy];
                }
            }
        }

        stationary_dist = new_dist;
    }
}



////////////////////////////////////////////////////// Block 6 //////////////////////////////////////////////////////
double aggregate_capital() {
    double total = 0.0;
    for (int iy = 0; iy < ny; ++iy) {
        for (int ia = 0; ia < na; ++ia) {
            total += stationary_dist[iy][ia] * agrid[policy_idx[iy][ia]];
        }
    }
    return total;
}


////////////////////////////////////////////////////// Block 6 //////////////////////////////////////////////////////
int main() {
    create_asset_grid();
    solve_value_function();
    compute_stationary_distribution();

    double K = aggregate_capital();
    cout << "Aggregate capital: " << K << endl;

    return 0;
}

