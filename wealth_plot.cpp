#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <cstdlib>
#include <string>

int main() {
    // Path to the CSV file
    std::string filePath = "output/wealth_distribution.csv";

    // Vectors to hold the asset values and densities
    std::vector<double> assets;
    std::vector<double> densities;

    // Variables to hold the input from the file
    double asset, density;

    // Open the file
    std::ifstream file(filePath);

    if (!file.is_open()) {
        std::cerr << "Failed to open file at " << filePath << std::endl;
        return 1;
    }

    // Skip the header line
    std::string header;
    getline(file, header);

    // Read the data from the file
    while (file >> asset >> density) {
        assets.push_back(asset);
        densities.push_back(density);

        // Skip the comma
        file.ignore(1);
    }

    file.close();

    // Display a simple histogram
    std::cout << "Histogram of Wealth Distribution:" << std::endl;
    for (size_t i = 0; i < densities.size(); i++) {
        std::cout << assets[i] << " | ";
        int barLength = static_cast<int>(densities[i] * 1000); // Scale factor for visualization
        for (int j = 0; j < barLength; j++) {
            std::cout << "*";
        }
        std::cout << std::endl;
    }

    return 0;
}

// g++ -std=c++11 wealth_plot.cpp -o wealth_plot
// ./wealth_plot
