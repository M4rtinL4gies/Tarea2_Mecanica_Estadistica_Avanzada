#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

// El parámetro nbr de molecular field, junto con los nbrs, la energía inicial y los índices para el vector de la energía deben ser cambiados cuando se cambien las condiciones de borde periodicas - no periodicas.

// Function declarations
void gray_flip(std::vector<int>& t, int N, int& k);
std::vector<std::array<int, 4>> nbrs_pbc(int N);
std::vector<std::vector<int>> nbrs_no_pbc(int N);
int molecular_field(const std::vector<int>& spins, const std::vector<std::array<int, 4>>& nbr, int k); // Change this for no_pbc
void write_densidades(int N, const std::vector<long long>& densities);
void write_densidades_m(int N, const std::vector<std::vector<long long>>& densities);

// Function definitions
void gray_flip(std::vector<int>& t, int N, int& k) {
    k = t[0];
    if (k > N) return;
    t[k - 1] = t[k];
    t[k] = k + 1;
    if (k != 1) t[0] = 1;
}

std::vector<long long> enumerate_ising(int N) {
    std::vector<long long> densities(4 * N + 1, 0);
    std::vector<int> spins(N, -1);
    std::vector<int> t(N + 1);
    for (int i = 0; i <= N; ++i) {
        t[i] = i + 1;
    }
    auto nbr = nbrs_pbc(N); // Change this for no_pbc
    int e = -2 * N;            // Change this to 2 * L * (L - 1), and change shift index to this value too
    densities[e + 2 * N] = 2;  // Shift index to fit in vector

    for (long long i = 0; i < (1ULL << (N - 1)) - 1; ++i) {
        int k;
        gray_flip(t, N, k);
        int h = molecular_field(spins, nbr, k - 1);

        e += 2 * spins[k - 1] * h;
        densities[e + 2 * N] += 2;

        spins[k - 1] *= -1;
    }

    return densities;
}

std::vector<std::vector<long long>> enumerate_ising_m(int N) {
    std::vector<std::vector<long long>> densities(4 * N + 1, std::vector<long long>(2 * N + 1, 0));  // Each row is an energy, each column is a magnetization
    std::vector<int> spins(N, -1);  // Initialize spins vector with -1
    std::vector<int> t(N + 1);
    for (int i = 0; i <= N; ++i) {
        t[i] = i + 1;
    }
    auto nbr = nbrs_pbc(N); // Change this for no_pbc
    int e = -2 * N; int m = - N;
    densities[e + 2 * N][m + N] = 1;  // Shift index to fit in matrix

    //std::cout << " " << std::endl;
    //std::cout << "Energy Value: " << e << ", #: " << densities[e + 2 * N] << std::endl;
    //std::cout << std::endl;

    for (long long i = 0; i < (1ULL << N) - 1; ++i) {
        //std::cout << "Iteration: " << i << std::endl;
        int k;
        gray_flip(t, N, k);
        int h = molecular_field(spins, nbr, k - 1);

        e += 2 * spins[k - 1] * h; m -= 2 * spins[k - 1];
        densities[e + 2 * N][m + N] += 1;
        //std::cout << "Energy Value: " << e << ", #: " << densities[e + 2 * N] << std::endl;

        spins[k - 1] *= -1;
        //std::cout << "spins: ";
        //for (const auto& element : spins){
        //    std::cout << element << " ";
        //}
        //std::cout << std::endl;
        //std::cout << std::endl;
    }

    return densities;
}

std::vector<std::array<int, 4>> nbrs_pbc(int N) {
    int L = std::sqrt(N);
    std::vector<std::array<int, 4>> nbr(N);
    for (int i = 0; i < N; ++i) {
        int right = ((i / L) * L + (i + 1) % L);
        int bottom = (i + L) % N;
        int left = ((i / L) * L + (i - 1 + L) % L);
        int top = (i - L + N) % N;
        nbr[i] = {right, bottom, left, top};
    }
    return nbr;
}

std::vector<std::vector<int>> nbrs_no_pbc(int N) {
    // In order to use this function, the nbr input in the molecular field function must be changed to a std::vector<std::vector<int>> type. Although this is a little bit inconvenient, the result is a more optimized program when using pbc
    int L = std::sqrt(N);  // Side length of the grid
    std::vector<std::vector<int>> nbr(N);

    for (int i = 0; i < N; ++i) {
        int row = i / L;  // Row index
        int col = i % L;  // Column index

        // Right
        if (col < L - 1) {
            nbr[i].push_back(i + 1);
        }
        // Bottom
        if (row < (N / L) - 1) {
            nbr[i].push_back(i + L);
        }
        // Left
        if (col > 0) {
            nbr[i].push_back(i - 1);
        }
        // Top
        if (row > 0) {
            nbr[i].push_back(i - L);
        }
    }
    return nbr;
}

// Change this for no_pbc
int molecular_field(const std::vector<int>& spins, const std::vector<std::array<int, 4>>& nbr, int k) {
    int sum = 0;
    const auto& nbrs_k = nbr[k];
    for (int n : nbrs_k) {
        sum += spins[n];
    }
    return sum;
}

void write_densidades(int N, const std::vector<long long>& densities) {
    // Change this for no pbc
    std::ofstream file("densidades_pbc.txt", std::ios::app);
    file << "E\t" << (int(std::sqrt(N))) << "x" << (int(std::sqrt(N))) << "\n";
    for (int i = -2 * N; i <= 2 * N; i += 4) {
        file << i << "\t" << densities[i + 2 * N] << "\n";
    }
    file << "\n";
}

void write_densidades_m(int N, const std::vector<std::vector<long long>>& densities) {
    std::ofstream file("densidades_em.txt", std::ios::app);
    file << (int(std::sqrt(N))) << "x" << (int(std::sqrt(N))) << "\n";
    file << "E\\M";
    for (int i = - N; i <= N; i += 2) {
        file << "\t" << i;
    }
    file << "\n";
    for (int i = -2 * N; i <= 2 * N; i += 4) {
        file << i;
        for (int j = - N; j <= N; j += 2) {
        file << "\t" << densities[i + 2 * N][j + N];
        }
        file << "\n";
    }
    file << "\n";
}

int main() {
    int N = 36;

    // Para las densidades únicamente energéticas
    //auto densities = enumerate_ising(N);
    //write_densidades(N, densities);

    // Para las densidades energéticas y magnéticas
    auto densities = enumerate_ising_m(N);
    write_densidades_m(N, densities);

    return 0;
}
