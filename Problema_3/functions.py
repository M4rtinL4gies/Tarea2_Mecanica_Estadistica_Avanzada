import numpy as np

# Functions
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

# Markov Ising
def markov_ising(spins: list, E, nbr, T):
    N = len(spins)
    rng = np.random.default_rng()

    k = rng.integers(0, N, endpoint = False)
    h = molecular_field(spins, nbr, k)
    delta_E = 2 * h * spins[k]
    gamma = np.exp(- delta_E / T)

    x = rng.uniform(0,1,1)[0]
    if x < gamma:
        spins[k] *= -1
        E += delta_E
    return spins, E

def calor_especifico(energy_list, T):
    E_mean = np.mean(energy_list)
    E2_mean = np.mean(np.square(energy_list))
    
    c_V = (E2_mean - E_mean**2) / T**2
    return c_V

def energy_ising(spins, nbr):
    E = 0
    for k, spin in enumerate(spins):
        nbrs_list = nbr[k]
        for i, n in enumerate(nbrs_list):
            E -= spin * spins[n]
    
    return int(E / 2)

def write_data(e_list, cv_list, T):
    with open("data/energy.txt", "a") as file:
        file.write("T\t<e>\tcv\n")
        for i in range(len(T)):
            file.write("{}\t{}\t{}\n".format(T[i], e_list[i], cv_list[i]))
        file.write("\n")

if __name__ == "__main__":
    N = 36
    L = int(np.sqrt(N))
    n_samples = 1e6
    t_equilibrio = 1000 # Nº de sweeps

    rng = np.random.default_rng()
    spins = rng.choice([-1,1], size = N)
    nbr = nbrs(N, True)
    E = energy_ising(spins, nbr)
    T = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]

    e_list_all_T = []
    m_list_all_T = []
    for t in T:
        e_list = []
        m_list = []
        for i in range(t_equilibrio * N): # Dejamos que pasen t_equilibrio sweeps
            spins, E = markov_ising(spins, E, nbr, t)

        # Ahora tomamos samples cada un sweep
        for i in range(n_samples):
            for j in range(N):
                spins, E = markov_ising(spins, E, nbr, t)
            e_list.append(E)
            m_list.append(sum(np.array(spins)))
        
        e_list_all_T.append(e_list)
        m_list_all_T.append(m_list)

    # Calculamos ahora la e media y el calor específico
    e_mean = []
    cv = []
    m_mean = []
    for i, list in enumerate(e_list_all_T):
        e_mean.append(np.mean(np.array(list))/N)
        cv.append(calor_especifico(np.array(list), T[i])/N)
        m_mean.append(np.mean(np.array(list))/N)

    # Guardamos los datos:
    write_data(e_mean, cv, T)