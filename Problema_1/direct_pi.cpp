#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <numeric>  // for accumulate
#include <string>   // for string
#include <sstream>  // for stringstream

const double PI = 3.14159265358979323846;

// Function to perform the Monte Carlo sampling to estimate Pi
double direct_pi(int N) {
    int n_hits = 0;
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);

    for (int i = 0; i < N; ++i) {
        double x = distribution(rng);
        double y = distribution(rng);
        if (x * x + y * y < 1) {
            ++n_hits;
        }
    }
    return static_cast<double>(n_hits) / N;
}

// Function to convert a double to a string
std::string to_string(double value) {
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

// Function to calculate mean, standard deviation, and mean squared deviation from Pi/4
std::vector<std::string> mean_and_stdv(const std::vector<double>& list) {
    double m = std::accumulate(list.begin(), list.end(), 0.0) / list.size();

    // Calculate standard deviation
    double stdv = 0.0;
    for (double value : list) {
        stdv += (value - m) * (value - m);
    }
    stdv = std::sqrt(stdv / list.size());

    // Calculate mean squared deviation from Pi/4
    double mcd = 0.0;
    for (double value : list) {
        mcd += (value - (PI / 4)) * (value - (PI / 4));
    }
    mcd /= list.size();

    // Convert the values to strings and return them
    return {to_string(m), to_string(stdv), to_string(mcd)};
}

// Function to generate data for different N and write it to a file
void make_direct_data(const std::vector<int>& runs, int n) {
    std::vector<std::vector<double>> results;
    std::vector<std::vector<std::string>> results_errors;

    for (int r : runs) {
        std::vector<double> results_n;
        for (int i = 0; i < n; ++i) {
            results_n.push_back(direct_pi(r));
        }
        results_errors.push_back(mean_and_stdv(results_n));
        results.push_back(results_n);
    }

    // Write results to the file
    std::ofstream file("data_direct_pi.txt");

    // Writing column headers
    file << "10**1\t10**2\t10**3\t10**4\t10**5\t10**6\t10**7\t10**8\n";
    
    // Writing individual results
    for (int i = 0; i < n; ++i) {
        for (size_t j = 0; j < runs.size(); ++j) {
            file << results[j][i] << "\t";
        }
        file << "\n";
    }
    
    file << "\n";

    // Writing mean, standard deviation, and mean squared deviation for each run
    for (size_t j = 0; j < runs.size(); ++j) {
        file << "1e" << (j + 1) << "\t" << results_errors[j][0] << "\t" 
             << results_errors[j][1] << "\t" << results_errors[j][2] << "\n";
    }

    file.close();
}

int main() {
    // Runs for different N values
    std::vector<int> runs = {10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000};
    int n = 20;  // Number of simulations for each N
    make_direct_data(runs, n);
    return 0;
}