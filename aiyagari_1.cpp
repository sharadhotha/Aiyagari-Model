#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

// Model parameters
const double beta = 0.96;   // Discount factor
const double alpha = 0.36;  // Capital share
const double delta = 0.08;  // Depreciation rate
const double sigma = 2.0;   // CRRA coefficient
const int grid_size = 100;  // Number of asset grid points
const int income_states = 2;// Number of income states
const double rho = 0.9;     // Persistence of income shocks
const double sigma_e = 0.1; // Std. dev. of shocks
const int max_iter = 1000;  // Maximum iterations
const double tol = 1e-6;    // Convergence tolerance

// Utility function (CRRA)
double utility(double c) {
    return (sigma == 1) ? log(c) : (pow(c, 1 - sigma) - 1) / (1 - sigma);
}

// Solve household problem using Value Function Iteration
void solve_value_function(vector<vector<double>> &V, vector<vector<int>> &policy, vector<double> &asset_grid, vector<double> &income_grid, double r, double w) {
    vector<vector<double>> V_new = V;
    double diff;
    int iter = 0;

    do {
        diff = 0.0;
        for (int i = 0; i < income_states; i++) {
            for (int j = 0; j < grid_size; j++) {
                double best_value = -1e9;
                int best_choice = 0;
                double a = asset_grid[j];
                double y = income_grid[i] * w;
                for (int k = 0; k < grid_size; k++) {
                    double a_next = asset_grid[k];
                    double c = (1 + r) * a + y - a_next;
                    if (c > 0) {
                        double val = utility(c) + beta * V[i][k];
                        if (val > best_value) {
                            best_value = val;
                            best_choice = k;
                        }
                    }
                }
                V_new[i][j] = best_value;
                policy[i][j] = best_choice;
                diff = max(diff, fabs(V_new[i][j] - V[i][j]));
            }
        }
        V = V_new;
        iter++;
    } while (diff > tol && iter < max_iter);
}

int main() {
    // Asset grid
    vector<double> asset_grid(grid_size);
    for (int i = 0; i < grid_size; i++) {
        asset_grid[i] = i * 0.1;
    }

    // Income states (simplified binary process)
    vector<double> income_grid = {0.5, 1.5};

    // Initialize value function and policy function
    vector<vector<double>> V(income_states, vector<double>(grid_size, 0.0));
    vector<vector<int>> policy(income_states, vector<int>(grid_size, 0));
    
    // Interest rate and wage from firm’s problem (assume fixed for now)
    double r = 0.04;
    double w = 1.0;
    
    // Solve the household’s problem
    solve_value_function(V, policy, asset_grid, income_grid, r, w);
    
    // Print policy function (optimal asset choice)
    cout << "Optimal savings policy (next period asset choice):" << endl;
    for (int i = 0; i < income_states; i++) {
        cout << "Income state " << i << ": ";
        for (int j = 0; j < grid_size; j += 10) { // Print every 10th asset point
            cout << asset_grid[policy[i][j]] << " ";
        }
        cout << endl;
    }
    
    return 0;
}
