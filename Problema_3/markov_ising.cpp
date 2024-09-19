#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <numeric>
#include <fstream>

// Function to initialize neighbors with periodic boundary conditions (PBC)
std::vector<std::vector<int>> nbrs(int N, bool pbc) {
    int L = static_cast<int>(std::sqrt(N));
    std::vector<std::vector<int>> nbr(N);

    if (pbc) {
        for (int i = 0; i < N; ++i) {
            nbr[i] = {
                (i / L) * L + (i + 1) % L,  // Right
                (i + L) % N,               // Bottom
                (i / L) * L + (i - 1 + L) % L, // Left
                (i - L + N) % N            // Top
            };
        }
    } else {
        for (int i = 0; i < N; ++i) {
            int row = i / L;
            int col = i % L;

            // Right
            if (col < L - 1) nbr[i].push_back(i + 1);
            // Bottom
            if (row < (N / L) - 1) nbr[i].push_back(i + L);
            // Left
            if (col > 0) nbr[i].push_back(i - 1);
            // Top
            if (row > 0) nbr[i].push_back(i - L);
        }
    }

    return nbr;
}

// Function to calculate the molecular field
int molecular_field(const std::vector<int>& spins, const std::vector<std::vector<int>>& nbr, int k) {
    int h = 0;
    for (int i : nbr[k]) {
        h += spins[i];
    }
    return h;
}

// Markov Ising function
void markov_ising(std::vector<int>& spins, int& E, const std::vector<std::vector<int>>& nbr, double T, std::mt19937& rng) {
    int N = spins.size();
    std::uniform_int_distribution<int> dist(0, N - 1);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    int k = dist(rng);
    int h = molecular_field(spins, nbr, k);
    int delta_E = 2 * h * spins[k];
    double gamma = std::exp(-delta_E / T);

    if (uniform(rng) < gamma) {
        spins[k] *= -1;
        E += delta_E;
    }
}

// Function to calculate specific heat
double calor_especifico(const std::vector<int>& energy_list, double T) {
    double E_mean = std::accumulate(energy_list.begin(), energy_list.end(), 0.0) / energy_list.size();
    double E2_mean = std::accumulate(energy_list.begin(), energy_list.end(), 0.0, [](double sum, int e) { return sum + e * e; }) / energy_list.size();
    return (E2_mean - E_mean * E_mean) / (T * T);
}

// Function to calculate energy of Ising configuration
int energy_ising(const std::vector<int>& spins, const std::vector<std::vector<int>>& nbr) {
    int E = 0;
    for (size_t k = 0; k < spins.size(); ++k) {
        for (int n : nbr[k]) {
            E -= spins[k] * spins[n];
        }
    }
    return E / 2; // To avoid double counting
}

// Function to write data
void write_data(const std::vector<double>& e_list, const std::vector<double>& cv_list, const std::vector<double>& T) {
    std::ofstream file("data/energy.txt", std::ios::app);
    file << "T\t<e>\tcv\n";
    for (size_t i = 0; i < T.size(); ++i) {
        file << T[i] << "\t" << e_list[i] << "\t" << cv_list[i] << "\n";
    }
    file << "\n";
    file.close();
}

// Function to write data magnetization
void write_data_m(const std::vector<double>& m_list, const std::vector<double>& T, int N) {
    std::ofstream file("data/magnetization.txt", std::ios::app);
    file << (int(std::sqrt(N))) << "x" << (int(std::sqrt(N))) << "\n";
    file << "T\t<|m|>\n";
    for (size_t i = 0; i < T.size(); ++i) {
        file << T[i] << "\t" << m_list[i] << "\n";
    }
    file << "\n";
    file.close();
}

int main() {
    int N = 1024;
    int L = static_cast<int>(std::sqrt(N));
    int n_samples = 1e5;
    int t_equilibrio = 1000;

    std::random_device rd;
    std::mt19937 rng(rd());

    std::vector<std::vector<int>> nbr = nbrs(N, true);
    std::vector<double> T;
    for (double t = 0; t <= 5.0; t += 0.1) {
        T.push_back(t);
    }

    std::vector<std::vector<int>> e_list_all_T(T.size());
    std::vector<std::vector<int>> m_list_all_T(T.size());

    for (size_t t_idx = 0; t_idx < T.size(); ++t_idx) {
        std::uniform_int_distribution<int> dist_spin(-1, 1);
        std::vector<int> spins(N);
        for (int& s : spins) {
            s = dist_spin(rng) == 0 ? 1 : -1; // Initialize spins to -1 or 1
        }
        int E = energy_ising(spins, nbr);
        
        double t = T[t_idx];
        for (int i = 0; i < t_equilibrio * N; ++i) { // Equilibration sweeps
            markov_ising(spins, E, nbr, t, rng);
        }

        for (int i = 0; i < n_samples; ++i) { // Sampling
            for (int j = 0; j < N; ++j) {
                markov_ising(spins, E, nbr, t, rng);
            }
            e_list_all_T[t_idx].push_back(E);
            int m = std::accumulate(spins.begin(), spins.end(), 0);
            m_list_all_T[t_idx].push_back(std::abs(m));
        }
    }

    // Calculate mean energy and specific heat
    std::vector<double> e_mean(T.size());
    std::vector<double> cv(T.size());
    std::vector<double> m_mean(T.size());
    for (size_t i = 0; i < T.size(); ++i) {
        e_mean[i] = std::accumulate(e_list_all_T[i].begin(), e_list_all_T[i].end(), 0.0) / N / e_list_all_T[i].size();
        cv[i] = calor_especifico(e_list_all_T[i], T[i]) / N;
        m_mean[i] = std::accumulate(m_list_all_T[i].begin(), m_list_all_T[i].end(), 0.0) / N / m_list_all_T[i].size();
    }

    // Write data
    write_data(e_mean, cv, T);
    write_data_m(m_mean, T, N);

    return 0;
}
