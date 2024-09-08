#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <numeric>
#include <string>
#include <sstream>

// Constants
const double PI = 3.14159265358979323846;

// Function to perform Markov Chain Monte Carlo sampling to estimate Pi
std::pair<double, double> markov_pi(int N, double delta) {
    int n_hits = 0;
    int n_rej = 0;
    double x = 1.0, y = 1.0;

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<double> dist(-delta, delta);

    for (int i = 0; i < N; ++i) {
        double d_x = dist(rng);
        double d_y = dist(rng);

        if (std::abs(x + d_x) < 1.0 && std::abs(y + d_y) < 1.0) {
            x += d_x;
            y += d_y;
        } else {
            ++n_rej;
        }

        if (x * x + y * y < 1.0) {
            ++n_hits;
        }
    }

    return std::make_pair(static_cast<double>(n_hits) / N, static_cast<double>(n_rej) / N);
}

// Function to convert double to string
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

// Function to generate Markov data for different N and write to file
void make_markov_data_n(const std::vector<int>& runs, int n) {
    std::vector<std::vector<double>> results;
    std::vector<std::vector<std::string>> results_errors;

    for (int r : runs) {
        std::vector<double> results_n;
        for (int i = 0; i < n; ++i) {
            results_n.push_back(markov_pi(r, 0.3).first);  // Taking the first element (hit ratio)
        }
        results_errors.push_back(mean_and_stdv(results_n));
        results.push_back(results_n);
    }

    // Write results to the file
    std::ofstream file("data_markov_pi.txt");

    // File header
    file << "For different N's...\n";
    file << "10**1\t10**2\t10**3\t10**4\t10**5\t10**6\t10**7\t10**8\n";

    // Write individual results
    for (int i = 0; i < n; ++i) {
        for (size_t j = 0; j < runs.size(); ++j) {
            file << results[j][i] << "\t";
        }
        file << "\n";
    }
    file << "\n";

    // Write mean, std deviation, and mean squared deviation for each N
    for (size_t j = 0; j < runs.size(); ++j) {
        file << "1e" << (j + 1) << "\t" << results_errors[j][0] << "\t" 
             << results_errors[j][1] << "\t" << results_errors[j][2] << "\n";
    }
    file << "\n";

    file.close();
}

// Function to generate Markov data for different deltas and write to file
void make_markov_data_delta(const std::vector<double>& deltas, int n, int N) {
    std::vector<std::vector<double>> results_delta;
    std::vector<std::vector<std::string>> results_delta_errors;
    std::vector<std::vector<double>> results_rej;
    std::vector<std::vector<std::string>> results_rej_errors;

    // Generate results for each delta
    for (double delta : deltas) {
        std::vector<double> results_d;
        std::vector<double> results_r;
        
        for (int i = 0; i < n; ++i) {
            auto [hits, rej] = markov_pi(N, delta);
            results_d.push_back(hits);
            results_r.push_back(rej);
        }

        results_delta_errors.push_back(mean_and_stdv(results_d));
        results_rej_errors.push_back(mean_and_stdv(results_r));
        results_delta.push_back(results_d);
        results_rej.push_back(results_r);
    }

    // Write results to the file
    std::ofstream file("data_markov_pi.txt", std::ios_base::app);  // Append mode

    // File header for different deltas
    file << "For different deltas...\n";
    file << "0.0 to 3.0 every 0.05\n";

    // Write individual results for different deltas
    for (int i = 0; i < n; ++i) {
        for (size_t j = 0; j < deltas.size(); ++j) {
            file << results_delta[j][i] << "\t";
        }
        file << "\n";
    }
    file << "\n";

    // Write mean, standard deviation, and mean squared deviation for each delta
    for (size_t j = 0; j < deltas.size(); ++j) {
        file << (j * 0.05) << "\t" << results_delta_errors[j][0] << "\t" 
             << results_delta_errors[j][1] << "\t" << results_delta_errors[j][2] << "\n";
    }
    file << "\n";

    // Write rejection rate
    file << "Tasa de rechazo...\n";
    for (int i = 0; i < n; ++i) {
        for (size_t j = 0; j < deltas.size(); ++j) {
            file << results_rej[j][i] << "\t";
        }
        file << "\n";
    }
    file << "\n";

    // Write rejection errors
    for (size_t j = 0; j < deltas.size(); ++j) {
        file << (j * 0.05) << "\t" << results_rej_errors[j][0] << "\t"
             << results_rej_errors[j][1] << "\t" << results_rej_errors[j][2] << "\n";
    }
    file << "\n";

    file.close();
}

int main() {
    // Runs for different N values
    std::vector<int> runs = {10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000};
    int n = 20;  // Number of simulations for each N
    //make_markov_data_n(runs, n);  // Original function call

    // Deltas for Markov data
    std::vector<double> deltas = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3.0};
    int N = 1000000;  // Number of steps for Markov chain
    make_markov_data_delta(deltas, n, N);  // Call for the new function

    return 0;
}

