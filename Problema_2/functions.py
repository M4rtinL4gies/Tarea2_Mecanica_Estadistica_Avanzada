# Imports
import numpy as np
from copy import deepcopy

# Functions
def all_configurations(k, t, N):
    conf_nueva = np.full(N, -1)
    configurations = np.array([conf_nueva])

    for _ in range(2**N-1):
        t, k = gray_flip(t, N)

        conf_nueva = configurations[-1].copy()
        print(conf_nueva)
        conf_nueva[k-1] = -1 * conf_nueva[k-1]
        configurations = np.append(configurations, [conf_nueva], axis = 0)

    return configurations

def write_configurations_file(configurations):
    with open("configurations.txt", "a") as file:
        file.write("i\t\t(s_1 ... s_{})".format(len(configurations[0])) + "\n")
        for i in range(len(configurations)):
            file.write(str(i+1) + "\t\t" + "\t".join("+" if j == 1 else "-" for j in configurations[i]) + "\n")
        file.write("\n")

def write_densidades(N, densidades):
    with open("densidades_pbc.txt", "a") as file:
        file.write("E\t{}x{}\n".format(int(N**(1/2)), int(N**(1/2))))
        for i in range(0, 2*N+1, 4):
            file.write("{}\t{}".format(i, densidades[i]) + "\n")
        file.write("\n")

def leer_estados():
    with open("data/densidades_em.txt", "r") as file:
        data = file.readlines()

    index_1 = data.index("\n")
    index_2 = len(data) - data[::-1][1:].index("\n") - 1

    data_2 = data[0:index_1]
    data_4 = data[index_1+1:index_2-1]
    data_6 = data[index_2:-1]
    datas = [data_2, data_4, data_6]

    dics = []
    for data in datas:
        mag = data[1].strip("\n").split("\t")[1:]
        N_E_M = {}
        # Process the rest of the lines to populate N(E, M)
        for line in data[2:]:
            values = list(map(int, line.strip().split()))
            E = values[0]  # First value in each line is the energy E
            N_E_M[E] = {int(M): count for M, count in zip(mag, values[1:])}
        dics.append(N_E_M)

    return dics[0], dics[1], dics[2]

def write_data_pi_M(data, T):
    with open("data/probabilities.txt", "a") as file:
        file.write("M\\T\t" + "\t".join([str(t) for t in T]) + "\n")
        for m in data:
            file.write("\t".join([str(p) for p in m]) + "\n")
        file.write("\n")

def write_data_binder(data, T, N):
    with open("data/binder.txt", "a") as file:
        file.write("T\tB(T)\tN = {}\n".format(N))
        for i in range(len(data)):
            file.write("{}\t{}\n".format(T[i], data[i]))
        file.write("\n")

def write_states(data):
    E = list(data.keys())
    M = list(data[0].keys())

    with open("data/densidades_em_2.txt", "a") as file:
        file.write("E\tM\tN\n")
        for e in E:
            for m in M:
                file.write("{}\t{}\t{}\n".format(e, m, data[e][m]))

        file.write("\n")

# ENUMERATE ISING
def gray_flip(t, N):
    k = t[0]
    if k > N: 
        return t, k
    t[k - 1] = t[k]
    t[k] = k + 1
    if k != 1: 
        t[0] = 1
    return t, k

def enumerate_ising(N):
    densidades = {e: 0 for e in range(-2*N, 2*N+1, 2)}
    spins = [-1 for i in range(N)]
    t = [i+1 for i in range(N+1)]
    e = -2*N
    densidades[e] = 2

    nbr = nbrs(N, True)

    for i in range(0, 2**(N-1) - 1):
        t, k = gray_flip(t, N)
        h = molecular_field(spins, nbr, k)
        e += 2*spins[k]*h
        
        try:
            densidades[e] += 2
        except KeyError as error:
            densidades[e] = 2
        
        spins[k] *= -1

    return densidades

def nbrs(N, pbc: bool):
    L = int(N**(1/2))
    if pbc:
        nbr = {
            i : ((i // L) * L + (i + 1) % L, (i + L) % N, (i // L) * L + (i - 1) % L, (i - L) % N) for i in range(N)
            } 
    else:
        nbr = {i: [] for i in range(N)}
        for i in range(N):
            row = i // L  # row
            col = i % L   # column

            # Right
            if col < L - 1:
                nbr[i].append(i + 1)
            # Bottom
            if row < (N // L) - 1:
                nbr[i].append(i + L)
            # Left
            if col > 0:
                nbr[i].append(i - 1)
            # Top
            if row > 0:
                nbr[i].append(i - L)

    return nbr
            
def molecular_field(spins, nbr, k):
    return sum([spins[i] for i in nbr[k]])


# PROBABILITIES AND BINDER
def calcular_pi_M(N_E_M, T):
    mag = list(N_E_M[0].keys())
    ener = list(N_E_M.keys())
    pi_M = [[m] for m in list(N_E_M[0].keys())]
    
    for t in T:
        pi_m = []
        beta = 1 / t
        Z = 0  # Función de partición

        for i, M in enumerate(mag):
            suma_E = 0
            for E in ener:
                suma_E += N_E_M[E][M] * np.exp(-beta * E)
            pi_m.append(suma_E)
            Z += suma_E

        for m,p in enumerate(pi_m):
            pi_M[m].append(p/Z)
        
    return pi_M

def binder(data_pi_m, T, N):
    m_vals = np.array([m[0]/N for m in data_pi_m])
    pi_vals = []
    for c in range(1,len(data_pi_m[0])):
        pi_m = [p[c] for p in data_pi_m]
        pi_vals.append(np.array(pi_m))
    
    m2 = []
    m4 = []

    for t in pi_vals:
        m2.append(np.sum(m_vals**2 * t))
        m4.append(np.sum(m_vals**4 * t))

    binder = []
    for i in range(len(m2)):
        binder.append((1/2) * (3 - m4[i] / m2[i]**2))

    return binder


if __name__ == "__main__":
    N = 36
    densidades = enumerate_ising(N)
    write_densidades(N, densidades)