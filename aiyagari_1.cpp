///////// Descrete set of possible Assets and possible Incomes

////////////////////////////////////////////////////// Block 1 //////////////////////////////////////////////////////
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

// CRRA Utility function
double utility(double c, double sigma) {
    if (c <= 0) return -1e10; // penalize infeasible consumption
    if (sigma == 1.0) return log(c);
    return pow(c, 1.0 - sigma) / (1.0 - sigma);
}


////////////////////////////////////////////////////// Block 2 //////////////////////////////////////////////////////
const int no_of_assets = 1000;                   // number of asset grid points
const int no_of_incomes = 2;                    // number of income states
const double beta = 0.96;                       // discount factor
const double sigma = 2.0;                       // risk aversion 
double r = 0.04;                                // interest rate
const double asset_min = 0.0;                   // borrowing constraint
const double asset_max = 5000.0;                  // max asset
const double tolerence = 1e-6;                  // tolerance
const int max_iterations = 1000;                // max iterations

const double alpha = 0.36;                      // capital share
const double depreciation = 0.08;               // depreciation

vector<double> assets_grid(no_of_assets);       // asset grid

void create_asset_grid() {
    for (int i = 0; i < no_of_assets; ++i) {
        assets_grid[i] = asset_min + i * (asset_max - asset_min) / (no_of_assets - 1);
    }
}


////////////////////////////////////////////////////// Block 3 //////////////////////////////////////////////////////
vector<double> incomes_set = {0.5, 1.5};
vector<vector<double> > transition_matrix = {{0.9, 0.1}, {0.1, 0.9}};





////////////////////////////////////////////////////// Block 4 //////////////////////////////////////////////////////
vector<vector<double> > Value_func(no_of_incomes, vector<double>(no_of_assets, 0.0));
vector<vector<int> > policy_idx(no_of_incomes, vector<int>(no_of_assets, 0));

void solve_value_function() {
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        double max_diff = 0.0;
        vector<vector<double> > Value_func_new = Value_func;

        for (int iy = 0; iy < no_of_incomes; ++iy) {
            for (int ia = 0; ia < no_of_assets; ++ia) {
                double max_val = -1e10;
                int best_a_idx = 0;

                for (int ja = 0; ja < no_of_assets; ++ja) {
                    double c = incomes_set[iy] + (1 + r) * assets_grid[ia] - assets_grid[ja];
                    double u = utility(c, sigma);

                    double EV = 0.0;
                    for (int jy = 0; jy < no_of_incomes; ++jy) {
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



////////////////////////////////////////////////////// Block 5 //////////////////////////////////////////////////////
vector<vector<double> > stationary_dist(no_of_incomes, vector<double>(no_of_assets, 0.0));

void compute_stationary_distribution() {
    // Initial guess: equal probability
    for (int iy = 0; iy < no_of_incomes; ++iy)
        for (int ia = 0; ia < no_of_assets; ++ia)
            stationary_dist[iy][ia] = 1.0 / (no_of_incomes * no_of_assets);

    for (int iter = 0; iter < 1000; ++iter) {
        vector<vector<double> > new_dist(no_of_incomes, vector<double>(no_of_assets, 0.0));

        for (int iy = 0; iy < no_of_incomes; ++iy) {
            for (int ia = 0; ia < no_of_assets; ++ia) {
                int ja = policy_idx[iy][ia];
                for (int jy = 0; jy < no_of_incomes; ++jy) {
                    new_dist[jy][ja] += stationary_dist[iy][ia] * transition_matrix[iy][jy];
                }
            }
        }

        stationary_dist = new_dist;
    }
}



////////////////////////////////////////////////////// Block 6 //////////////////////////////////////////////////////
double aggregate_capital() {
    double total = 0.0;
    for (int iy = 0; iy < no_of_incomes; ++iy) {
        for (int ia = 0; ia < no_of_assets; ++ia) {
            total += stationary_dist[iy][ia] * assets_grid[policy_idx[iy][ia]];
        }
    }
    return total;
}

////////////////////////////////////////////////////// Block 8 //////////////////////////////////////////////////////
double household_aggregate_capital(double r) {
    // Update global r
    ::r = r;

    solve_value_function();
    compute_stationary_distribution();
    return aggregate_capital();
}

////////////////////////////////////////////////////// Block 9 ////////////////////////////////////////////////////// Solves for r
void solve_general_equilibrium() {
    double r_low = 0.005;
    double r_high = 0.04;
    double tol_r = 1e-4;
    double K_supply, K_demand;

    while (r_high - r_low > tol_r) {
        double r_mid = 0.5 * (r_low + r_high);
        r = r_mid;

        K_supply = household_aggregate_capital(r_mid);
        K_demand = pow((alpha / (r_mid + depreciation)), 1.0 / (1.0 - alpha));

        if (K_supply > K_demand) {
            r_high = r_mid;
        } else {
            r_low = r_mid;
        }
    }

    r = 0.5 * (r_low + r_high);
    K_supply = household_aggregate_capital(r);

    cout << "Equilibrium interest rate r: " << r << endl;
    cout << "Equilibrium capital: " << K_supply << endl;
}

////////////////////////////////////////////////////// Block 10 ////////////////////////////////////////////////////// Capital vs r
void generate_supply_demand_data(const string& filename) {
    ofstream file(filename);
    file << "r,K_supply,K_demand\n";

    for (double r_test = 0.005; r_test <= 0.04; r_test += 0.001) {
        double K_supply = household_aggregate_capital(r_test);
        double K_demand = pow(alpha / (r_test + depreciation), 1.0 / (1.0 - alpha));
        file << r_test << "," << K_supply << "," << K_demand << "\n";
    }

    file.close();
    cout << "Saved supply/demand data to: " << filename << endl;
}

void generatePlot(const std::string& csvFileName) {
    // Create the subfolder if it doesn't exist
    system("mkdir -p plots"); // "mkdir -p" ensures the folder is created

    // Create the Gnuplot command string
    std::string gnuplotCommand = "gnuplot -e \"set terminal png size 800,600; "
                                 "set output 'plots/plot.png'; "
                                 "set title 'CSV Data Plot'; "
                                 "set xlabel 'Interest Rate (r)'; "
                                 "set ylabel 'Capital (K)'; "
                                 "set datafile separator ','; "
                                 "plot '" + csvFileName + "' using 1:2 with linespoints title 'K_{Supply}', "
                                 "'" + csvFileName + "' using 1:3 with linespoints title 'K_{Demand}'\"";

    // Execute the Gnuplot command
    int result = system(gnuplotCommand.c_str());
    if (result != 0) {
        std::cerr << "Error executing Gnuplot!" << std::endl;
        return;
    }

    std::cout << "Plot saved as 'plots/r_K_plot.png'" << std::endl;
}



////////////////////////////////////////////////////// Block 11 ////////////////////////////////////////////////////// Wealth Dist

void export_wealth_distribution(const string& filename) {
    ofstream file(filename);
    file << "asset,density\n";

    for (int ia = 0; ia < no_of_assets; ++ia) {
        double density = 0.0;
        for (int iy = 0; iy < no_of_incomes; ++iy) {
            density += stationary_dist[iy][ia];
        }
        file << fixed << setprecision(6) << assets_grid[ia] << "," << density << "\n";
    }

    file.close();
    cout << "Wealth distribution saved to: " << filename << endl;
}

void generateDensityPlot(const std::string& csvFileName) {
    // Create the subfolder if it doesn't exist
    system("mkdir -p plots"); // "mkdir -p" ensures the folder is created

    // Create the Gnuplot command string for plotting a density function
    std::string gnuplotCommand = "gnuplot -e \"set terminal png size 800,600; "
                                 "set output 'plots/density_plot.png'; "
                                 "set title 'Wealth Distribution'; "
                                 "set xlabel 'Wealth'; "
                                 "set ylabel 'Density'; "
                                 "set datafile separator ','; "
                                 "plot '" + csvFileName + "' using 1:2 with lines title 'Density'\"";

    // Execute the Gnuplot command
    int result = system(gnuplotCommand.c_str());
    if (result != 0) {
        std::cerr << "Error executing Gnuplot!" << std::endl;
        return;
    }

    std::cout << "Density plot saved as 'plots/wealth_dist.png'" << std::endl;
}



////////////////////////////////////////////////////// Block 12 //////////////////////////////////////////////////////
int main() {
    create_asset_grid();
    solve_general_equilibrium();
    generate_supply_demand_data("output/supply_demand.csv");
    generatePlot("output/supply_demand.csv");
    export_wealth_distribution("output/wealth_distribution.csv");
    generateDensityPlot("output/wealth_distribution.csv");

    return 0;
}


// g++ -std=c++11 aiyagari_1.cpp -o aiyagari_1
// ./aiyagari_1
// g++ -std=c++11 aiyagari_1.cpp -o aiyagari_1 \
  -I/opt/homebrew/include -L/opt/homebrew/lib \
  -lboost_iostreams -lboost_system

