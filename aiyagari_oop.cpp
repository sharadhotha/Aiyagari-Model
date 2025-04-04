#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <cstdlib>
#include <string>

using namespace std;

// Global Variables
double r = 0.04; // interest rate

class Household {
private:
    vector<double> assets_grid;
    vector<double> incomes;
    vector<vector<double>> transition_matrix;
    vector<vector<double>> Value_func;
    vector<vector<int>> policy_idx;
    vector<vector<double>> stationary_dist;
    double beta;
    double sigma;

public:
    Household(int no_of_assets, double asset_min, double asset_max, const vector<double>& incomes_set, const vector<vector<double>>& transition_matrix, double beta, double sigma) :
        incomes(incomes_set), transition_matrix(transition_matrix), beta(beta), sigma(sigma) {
        assets_grid.resize(no_of_assets);
        for (int i = 0; i < no_of_assets; ++i) {
            assets_grid[i] = asset_min + i * (asset_max - asset_min) / (no_of_assets - 1);
        }
        Value_func.resize(incomes.size(), vector<double>(no_of_assets, 0.0));
        policy_idx.resize(incomes.size(), vector<int>(no_of_assets, 0));
        stationary_dist.resize(incomes.size(), vector<double>(no_of_assets, 0.0));
    }

    double utility(double c) {
        if (c <= 0) return -1e10; // penalize infeasible consumption
        if (sigma == 1.0) return log(c);
        return pow(c, 1.0 - sigma) / (1.0 - sigma);
    }

    void solve_value_function(double tolerence, int max_iterations) {
        for (int iteration = 0; iteration < max_iterations; ++iteration) {
            double max_diff = 0.0;
            vector<vector<double>> Value_func_new = Value_func;

            for (size_t iy = 0; iy < incomes.size(); ++iy) {
                for (size_t ia = 0; ia < assets_grid.size(); ++ia) {
                    double max_val = -1e10;
                    int best_a_idx = 0;

                    for (size_t ja = 0; ja < assets_grid.size(); ++ja) {
                        double c = incomes[iy] + (1 + r) * assets_grid[ia] - assets_grid[ja];
                        double u = utility(c);
                        double EV = 0.0;

                        for (size_t jy = 0; jy < incomes.size(); ++jy) {
                            EV += transition_matrix[iy][jy] * Value_func[jy][ja];
                        }

                        double val = u + beta * EV;
                        if (val > max_val) {
                            max_val = val;
                            best_a_idx = ja;
                        }
                    }

                    Value_func_new[iy][ia] = max_val;
                    policy_idx[iy][ia] = best_a_idx;
                    max_diff = max(max_diff, fabs(Value_func_new[iy][ia] - Value_func[iy][ia]));
                }
            }

            Value_func = Value_func_new;
            if (max_diff < tolerence) break;
        }
    }

    void compute_stationary_distribution(int max_iterations) {
        // Initial guess: equal probability
        for (size_t iy = 0; iy < incomes.size(); ++iy)
            for (size_t ia = 0; ia < assets_grid.size(); ++ia)
                stationary_dist[iy][ia] = 1.0 / (incomes.size() * assets_grid.size());

        for (int iter = 0; iter < max_iterations; ++iter) {
            vector<vector<double>> new_dist(incomes.size(), vector<double>(assets_grid.size(), 0.0));

            for (size_t iy = 0; iy < incomes.size(); ++iy) {
                for (size_t ia = 0; ia < assets_grid.size(); ++ia) {
                    int ja = policy_idx[iy][ia];
                    for (size_t jy = 0; jy < incomes.size(); ++jy) {
                        new_dist[jy][ja] += stationary_dist[iy][ia] * transition_matrix[iy][jy];
                    }
                }
            }

            stationary_dist = new_dist;
        }
    }

    double aggregate_capital() {
        double total = 0.0;
        for (size_t iy = 0; iy < incomes.size(); ++iy) {
            for (size_t ia = 0; ia < assets_grid.size(); ++ia) {
                total += stationary_dist[iy][ia] * assets_grid[policy_idx[iy][ia]];
            }
        }
        return total;
    }
};

class Firm {
private:
    double alpha;
    double depreciation;

public:
    Firm(double alpha, double depreciation) : alpha(alpha), depreciation(depreciation) {}

    double capital_demand(double interest_rate) {
        return pow(alpha / (interest_rate + depreciation), 1.0 / (1.0 - alpha));
    }
};

void solve_general_equilibrium(Household& household, Firm& firm, double r_low, double r_high, double tol_r) {
    double K_supply, K_demand;

    while (r_high - r_low > tol_r) {
        double r_mid = 0.5 * (r_low + r_high);
        r = r_mid;

        K_supply = household.aggregate_capital();
        K_demand = firm.capital_demand(r_mid);

        if (K_supply > K_demand) {
            r_high = r_mid;
        } else {
            r_low = r_mid;
        }
    }

    r = 0.5 * (r_low + r_high);
    K_supply = household.aggregate_capital();

    cout << "Equilibrium interest rate r: " << r << endl;
    cout << "Equilibrium capital: " << K_supply << endl;
}

int main() {
    const int no_of_assets = 50000;
    const double asset_min = 0.0;
    const double asset_max = 5000.0;
    vector<double> incomes_set = {0.5, 1.5};
    vector<vector<double>> transition_matrix = {{0.9, 0.1}, {0.1, 0.9}};
    double beta = 0.96;
    double sigma = 2.0;
    double tolerence = 1e-6;
    int max_iterations = 1000;
    double alpha = 0.36;
    double depreciation = 0.08;

    Household household(no_of_assets, asset_min, asset_max, incomes_set, transition_matrix, beta, sigma);
    Firm firm(alpha, depreciation);

    household.solve_value_function(tolerence, max_iterations);
    household.compute_stationary_distribution(max_iterations);

    solve_general_equilibrium(household, firm, 0.005, 0.04, 1e-4);

    return 0;
}
