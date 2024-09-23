#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <numeric>
#include <fstream>
#include <omp.h>  // Es necesario incluir OpenMP

// Función para inicializar los vecinos con condiciones de frontera periódicas (PBC)
std::vector<std::vector<int>> nbrs(int N, bool pbc) {
    int L = static_cast<int>(std::sqrt(N));
    std::vector<std::vector<int>> nbr(N);

    if (pbc) {
        for (int i = 0; i < N; ++i) {
            nbr[i] = {
                (i / L) * L + (i + 1) % L,  // Derecha
                (i + L) % N,                // Abajo
                (i / L) * L + (i - 1 + L) % L, // Izquierda
                (i - L + N) % N             // Arriba
            };
        }
    } else {
        for (int i = 0; i < N; ++i) {
            int row = i / L;
            int col = i % L;

            // Derecha
            if (col < L - 1) nbr[i].push_back(i + 1);
            // Abajo
            if (row < (N / L) - 1) nbr[i].push_back(i + L);
            // Izquierda
            if (col > 0) nbr[i].push_back(i - 1);
            // Arriba
            if (row > 0) nbr[i].push_back(i - L);
        }
    }

    return nbr;
}

// Función para calcular el campo molecular
int molecular_field(const std::vector<int>& spins, const std::vector<std::vector<int>>& nbr, int k) {
    int h = 0;
    for (int i : nbr[k]) {
        h += spins[i];
    }
    return h;
}

// Función de Markov para el modelo Ising
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

// Función para calcular el calor específico
double calor_especifico(const std::vector<int>& energy_list, double T) {
    double E_mean = std::accumulate(energy_list.begin(), energy_list.end(), 0.0) / energy_list.size();
    double E2_mean = std::accumulate(energy_list.begin(), energy_list.end(), 0.0, [](double sum, int e) { return sum + e * e; }) / energy_list.size();
    return (E2_mean - E_mean * E_mean) / (T * T);
}

// Función para calcular la energía del sistema de Ising
int energy_ising(const std::vector<int>& spins, const std::vector<std::vector<int>>& nbr) {
    int E = 0;
    for (size_t k = 0; k < spins.size(); ++k) {
        for (int n : nbr[k]) {
            E -= spins[k] * spins[n];
        }
    }
    return E / 2; // Evitar doble conteo
}

// Función para escribir los datos de energía y calor específico
void write_data(const std::vector<double>& e_list, const std::vector<double>& cv_list, const std::vector<double>& T) {
    std::ofstream file("data/energy.txt", std::ios::app);
    file << "T\t<e>\tcv\n";
    for (size_t i = 0; i < T.size(); ++i) {
        file << T[i] << "\t" << e_list[i] << "\t" << cv_list[i] << "\n";
    }
    file << "\n";
    file.close();
}

// Función para escribir los datos de magnetización
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

// Función para correr la simulación
void run_simulation(double temp, int N, int n_samples, int t_equilibrio, const std::vector<std::vector<int>>& nbr, std::vector<double>& e_mean, std::vector<double>& cv, std::vector<double>& m_mean, int idx) {
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> dist_spin(-1, 1);
    
    std::vector<int> spins(N);
    for (int& s : spins) {
        s = dist_spin(rng) == 0 ? 1 : -1; // Inicializar spins a -1 o 1
    }

    int E = energy_ising(spins, nbr);

    // Equilibración
    for (int i = 0; i < t_equilibrio * N; ++i) {
        markov_ising(spins, E, nbr, temp, rng);
    }

    std::vector<int> e_list, m_list;
    for (int i = 0; i < n_samples; ++i) {
        for (int j = 0; j < N; ++j) {
            markov_ising(spins, E, nbr, temp, rng);
        }
        // e_list.push_back(E);
        int m = std::accumulate(spins.begin(), spins.end(), 0);
        m_list.push_back(std::abs(m));
    }

    // e_mean[idx] = std::accumulate(e_list.begin(), e_list.end(), 0.0) / N / e_list.size();
    // cv[idx] = calor_especifico(e_list, temp) / N;
    m_mean[idx] = std::accumulate(m_list.begin(), m_list.end(), 0.0) / N / m_list.size();
}

int main() {
    int N = 1024;
    int n_samples = 1e5;
    int t_equilibrio = 1000;

    std::vector<std::vector<int>> nbr = nbrs(N, true);
    std::vector<double> T;
    for (double t = 0.1; t <= 5.0; t += 0.075) {
        T.push_back(t);
    }

    std::vector<double> e_mean(T.size()), cv(T.size()), m_mean(T.size());

    // Paralelización con OpenMP
    #pragma omp parallel for
    for (int t_idx = 0; t_idx < T.size(); ++t_idx) {
        run_simulation(T[t_idx], N, n_samples, t_equilibrio, nbr, e_mean, cv, m_mean, t_idx);
    }

    // Escribir los resultados
    // write_data(e_mean, cv, T);
    write_data_m(m_mean, T, N);

    return 0;
}
