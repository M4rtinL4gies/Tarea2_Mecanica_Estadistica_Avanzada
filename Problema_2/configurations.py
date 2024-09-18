"""
MÓDULO PARA OBTENER LAS POSIBLES CONFIGURACIONES DE N SPINS MEDIANTE GRAY-FLIP
"""

from functions import all_configurations, write_configurations_file
import numpy as np
import time

N = 4
t = np.array([i for i in range(1, N+2)])
start = time.time()
conf_table = all_configurations(t[0], t, N)
write_configurations_file(conf_table)

end = time.time()
print("Elapsed (with compilation) = %s" % (end - start))