import numpy as np

def direct_pi(N):
    """
    Sampling directo Montecarlo
    """
    n_hits = 0
    rng = np.random.default_rng()
    for _ in range(N):
        x = rng.uniform(-1,1,1)
        y = rng.uniform(-1,1,1)
        if (x**2 + y**2 < 1):
            n_hits += 1
    return float(n_hits / N)

def markov_pi(N, delta):
    """
    Sampling Markov Chain
    """
    n_hits = 0
    n_rej = 0
    x, y = 1, 1
    rng = np.random.default_rng()
    for _ in range(N):
        d_x = rng.uniform(-delta, delta, 1)
        d_y = rng.uniform(-delta, delta, 1)
        if abs(x + d_x) < 1 and abs(y + d_y) < 1:
            x += d_x
            y += d_y
        else:
            n_rej += 1
        if x**2 + y**2 < 1:
            n_hits += 1
    return float(n_hits/N), float(n_rej/N)

def mean_and_stdv(list):
    """
    Retorna parámetros importantes respecto a cada medición
    """
    m = sum(list)/len(list)
    stdv = np.std(np.array(list), dtype=np.float64)

    mcd = 0
    for i in np.array(list):
        mcd += (i - np.pi / 4)**2
    mcd /= len(list)

    return str(m), str(stdv), str(mcd)

def make_direct_data(runs, n):
    """
    Crea un archivo con len(runs) columnas, donde cada columna representa un N distinto. Los datos contienen n mediciones cada uno.
    El archivo incluye además datos de la media y desviación para cada N.
    """
    results = []
    results_errors = []

    for r in runs:
        results_n = np.array([])
        for i in range(n):
            results_n = np.append(results_n, direct_pi(r))

        results_errors.append(mean_and_stdv(results_n))
        results.append(results_n)


    with open("direct_pi.txt", "w") as file:
        file.write("10**1\t10**2\t10**3\t10**4\t10**5\t10**6\t10**7\t10**8\n")
        for i in range(n):
            for j in range(len(runs)):
                file.write(str(results[j][i]) + "\t")
            file.write("\n")
        file.write("\n")

        for j in range(len(runs)):
            file.write("10e{}\t".format(j+1) + "\t".join(results_errors[j]) + "\n")

def make_markov_data_n(runs, n):
    """
    Crea un archivo con len(runs) columnas, donde cada columna representa un N distinto. Los datos contienen n mediciones cada uno.
    El archivo incluye además datos de la media y desviación para cada N.
    """
    results = []
    results_errors = []

    for r in runs:
        results_n = np.array([])
        for i in range(n):
            results_n = np.append(results_n, markov_pi(r, 0.3)[0])

        results_errors.append(mean_and_stdv(results_n))
        results.append(results_n)


    with open("data_markov_pi.txt", "w") as file:
        file.write("For differents N's...\n")
        file.write("10**1\t10**2\t10**3\t10**4\t10**5\t10**6\t10**7\t10**8\n")
        for i in range(n):
            for j in range(len(runs)):
                file.write(str(results[j][i]) + "\t")
            file.write("\n")
        file.write("\n")

        for j in range(len(runs)):
            file.write("1e{}\t".format(j+1) + "\t".join(results_errors[j]) + "\n")
        file.write("\n")

def make_markov_data_delta(runs, n, N):
    """
    Crea un archivo con len(runs) columnas, donde cada columna representa un delta distinto. Los datos contienen n mediciones cada uno.
    El archivo incluye además datos de la media y desviación para cada delta.
    """
    results_delta = []
    results_delta_errors = []
    results_rej = []
    results_rej_errors = []

    for r in runs:
        results_d = np.array([])
        results_r = np.array([])
        for i in range(n):
            hits, rej = markov_pi(N, r)
            results_d = np.append(results_d, hits)
            results_r = np.append(results_r, rej)

        results_delta_errors.append(mean_and_stdv(results_d))
        results_rej_errors.append(mean_and_stdv(results_r))
        results_delta.append(results_d)
        results_rej.append(results_r)

    with open("data_markov_pi.txt", "a") as file:
        file.write("For differents deltas...\n")
        file.write("0.0\t0.2\t0.4\t0.6\t0.8\t1.0\t1.2\t1.4\t1.6\t1.8\t2.0\t2.2\t2.4\t2.6\t2.8\t3.0\n")
        for i in range(n):
            for j in range(len(runs)):
                file.write(str(results_delta[j][i]) + "\t")
            file.write("\n")
        file.write("\n")

        for j in range(len(runs)):
            file.write("{}\t".format(j * 0.2) + "\t".join(results_delta_errors[j]) + "\n")
        file.write("\n")

        file.write("Tasa de rechazo...\n")
        for i in range(n):
            for j in range(len(runs)):
                file.write(str(results_rej[j][i]) + "\t")
            file.write("\n")
        file.write("\n")

        for j in range(len(runs)):
            file.write("{}\t".format(j * 0.2) + "\t".join(results_rej_errors[j]) + "\n")
        file.write("\n")
