#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <queue>
#include <numeric>
#include <algorithm>

// Definir las constantes del sistema
const int L = 64;  // Tamaño de la red
const int N = L * L;  // Número total de spins

// Función para generar la lista de vecinos con condiciones de borde periódicas
std::vector<std::vector<int>> nbrs(int N, int L) {
    std::vector<std::vector<int>> nbr(N);
    for (int i = 0; i < N; ++i) {
        nbr[i] = {
            (i / L) * L + (i + 1) % L,  // Right
            (i + L) % N,               // Bottom
            (i / L) * L + (i - 1 + L) % L, // Left
            (i - L + N) % N            // Top
        };
    }
    return nbr;
}

// Función para calcular la energía del sistema
int energy_ising(const std::vector<int>& spins, const std::vector<std::vector<int>>& nbr) {
    int E = 0;
    for (size_t k = 0; k < spins.size(); ++k) {
        for (int n : nbr[k]) {
            E -= spins[k] * spins[n];
        }
    }
    return E / 2; // To avoid double counting
}

// Función para calcular el cambio de energía por el flip de un spin
int molecular_field(const std::vector<int>& spins, const std::vector<std::vector<int>>& nbr, int k) {
    int h = 0;
    for (int i : nbr[k]) {
        h += spins[i];
    }
    return h;
}

// Algoritmo de Wolff para generar clusters
void cluster_ising(std::vector<int>& spins, int N, const std::vector<std::vector<int>>& nbrs, double p_add, std::mt19937& rng, double& E, int& M) {
    std::uniform_int_distribution<int> randint(0, N - 1);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    int j = randint(rng);
    std::queue<int> pocket;
    pocket.push(j);

    std::vector<bool> visited(N, false);
    visited[j] = true;

    int initial_spin = spins[j];
    std::vector<int> cluster;  // Vector to store indices of spins to be flipped
    cluster.push_back(j);       // Add the initial spin index to the cluster

    // Cluster generation
    while (!pocket.empty()) {
        int k = pocket.front();
        pocket.pop();

        for (int l : nbrs[k]) {
            if (!visited[l] && spins[l] == initial_spin && uniform(rng) < p_add) {
                pocket.push(l);
                visited[l] = true;
                cluster.push_back(l);  // Add to the cluster
            }
        }
    }

    // Flip all spins in the cluster and calculate energy change
    for (int index : cluster) {
        spins[index] *= -1;  // Flip the spin
        // int h = molecular_field(spins, nbrs, index);  // Calculate molecular field for the flipped spin
        // E -= 2 * h * spins[index];  // Update energy for the flipped spin
        M += 2 * spins[index];
    }
}


// Funciones para escribir los resultados en un archivo
void write_data(const std::vector<double>& T_list, const std::vector<double>& e_mean_list, const std::vector<double>& cv_list) {
    std::ofstream file("data/energy.txt", std::ios::app);
    if (file.is_open()) {
        file << "T\t<E/N>\tCv/N\n";
        for (size_t i = 0; i < T_list.size(); ++i) {
            file << T_list[i] << "\t" << e_mean_list[i] << "\t" << cv_list[i] << "\n";
        }
        file << "\n";
        file.close();
        std::cout << "Data successfully written to 'data/cluster.txt'.\n";
    } else {
        std::cerr << "Error opening file for writing!\n";
    }
}

void write_data_m(const std::vector<double>& m_list, const std::vector<double>& b_list, const std::vector<double>& T, int N) {
    std::ofstream file("data/magnetization.txt", std::ios::app);
    file << (int(std::sqrt(N))) << "x" << (int(std::sqrt(N))) << "\n";
    file << "T\t<|m|>\tB(T)\n";
    for (size_t i = 0; i < T.size(); ++i) {
        file << T[i] << "\t" << m_list[i] << "\t" << b_list[i] << "\n";
    }
    file << "\n";
    file.close();
}

// Función para realizar las mediciones
void simulate_wolff(int L, int N, int n_samples, const std::vector<double>& T_list) {
    std::random_device rd;
    std::mt19937 rng(rd());  
    std::uniform_int_distribution<int> binary(0, 1);

    std::vector<std::vector<int>> neighbors = nbrs(N, L);

    std::vector<double> e_mean_list;
    std::vector<double> cv_list;
    std::vector<double> m_mean_list;
    std::vector<double> binder_list;

    for (double T : T_list) {
        double beta = 1.0 / T;
        double p_add = 1 - exp(-2.0 * beta);
        std::vector<int> spins(N);
        for (int k = 0; k < N; ++k) {
            spins[k] = binary(rng) == 0 ? 1 : -1;
        }

        double E = energy_ising(spins, neighbors);
        int M = std::accumulate(spins.begin(), spins.end(), 0);

        double Et = E;
        double Et2 = E * E;
        double Mt = std::abs(M);
        double Mt2 = M * M;
        double Mt4 = M * M * M * M;
        double m = static_cast<double>(M) / N;  // Calculate magnetization
        double mt = std::abs(m);                       // Sum |m|
        double mt2 = m * m;                            // Sum m^2
        double mt4 = m * m * m * m;                    // Sum m^4

        for (int i = 0; i < 1000; ++i) {
            cluster_ising(spins, N, neighbors, p_add, rng, E, M);
        }
        // Inside your main simulation loop
        for (int i = 0; i < n_samples; ++i) {
            cluster_ising(spins, N, neighbors, p_add, rng, E, M);
            
            // Accumulate observables
            double m = static_cast<double>(M) / N;  // Calculate magnetization
            mt += std::abs(m);                       // Sum |m|
            mt2 += m * m;                            // Sum m^2
            mt4 += m * m * m * m;                    // Sum m^4
        }

        // Now compute averages
        double m_mean = mt / (n_samples + 1);  // Average |M|, normalized by N
        double m2_mean = mt2 / (n_samples + 1);    // Average M^2
        double m4_mean = mt4 / (n_samples + 1);    // Average M^4

        // Correct Binder cumulant calculation:
        double binder = (1.0 / 2.0) * (3.0 - m4_mean / (m2_mean * m2_mean));

        // Store results
        m_mean_list.push_back(m_mean);
        binder_list.push_back(binder);

        std::cout << "Progreso " << T << "/" << T_list.back() << std::endl;
    }

    // write_data(T_list, e_mean_list, cv_list);
    write_data_m(m_mean_list, binder_list, T_list, N);
}

int main() {
    int n_samples = 1e5;
    std::vector<double> T_list;

    for (double T = 2.15; T <= 2.35; T += 0.025) {
        T_list.push_back(T);
    }

    simulate_wolff(L, N, n_samples, T_list);

    return 0;
}