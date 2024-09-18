"""
MÓDULO PARA OBTENER LA DISTRIBUCIÓN DE PROBABILIDADES DE LA MAGNETIZACIÓN
"""
# Imports
import numpy as np
from functions import leer_estados, calcular_pi_M, binder, write_data_pi_M, write_data_binder, write_states

if __name__ == "__main__":
    states_2, states_4, states_6 = leer_estados()

    T1 = [2.5, 3, 5]
    pi_M = calcular_pi_M(states_6, T1)
    write_data_pi_M(pi_M, T1)
    
    T2 = np.linspace(0.2, 50, 500)
    pi_M_2, pi_M_4, pi_M_6 = calcular_pi_M(states_2, T2), calcular_pi_M(states_4, T2), calcular_pi_M(states_6, T2)
    binder_2, binder_4, binder_6 = binder(pi_M_2, T2, 4), binder(pi_M_4, T2, 16), binder(pi_M_6, T2, 36)
    write_data_binder(binder_6, T2, 36)