import scqubits as scq
import qutip as qt
import numpy as np
import scipy as sp
from tqdm import tqdm

def get_transitions(qubit, evals_count):
    """
    Calculate the energy transitions for a given qubit.

    Parameters:
    qubit (Qubit object): The qubit for which transitions are to be calculated.
    evals_count (int): The number of eigenvalues to consider.

    Returns:
    ndarray: A 2D array of energy transitions between different eigenstates.
    """
    # Get the eigenvalues of the qubit
    eigenvals = qubit.eigenvals(evals_count)
    # Calculate the differences between each pair of eigenvalues
    return np.subtract.outer(eigenvals, eigenvals)

def compute_spectrum(system, fluxonium, phi_ext_values = np.linspace(-0.6,0.6,100), max_levels = 10):
    """
    Compute the transition energies for a range of external flux values.

    Parameters:
    - phi_ext_values: An array of external flux values to iterate over.
    - max_levels: The maximum energy level to consider for transitions.
    - system: The quantum system for which to calculate transitions.
    - get_transitions: A function that calculates transition energies given the system.

    Returns:
    - Dictionary containing transition energies for each pair of levels.
    """
    energy_level_transitions = {(i, j): [] for i in range(max_levels + 1) for j in range(i + 1, max_levels + 1)}

    for phi in tqdm(phi_ext_values):
        fluxonium.flux = phi  # Assuming the system has a 'flux' attribute
        transition_energies = get_transitions(system, max_levels + 1)

        for (i, j) in energy_level_transitions.keys():
            energy_diff = transition_energies[j, i]
            energy_level_transitions[(i, j)].append(energy_diff)

    return energy_level_transitions
def compute_chi_terms(fluxonium, ChainMode, max_level ):
    """
    Compute the terms in the square brackets for the chi quantity of a fluxonium system.

    Parameters:
    - fluxonium: A Fluxonium object from scqubits.
    - ChainMode: A chain mode custom object
    - max_level: The maximum level 'l' to consider in the fluxonium system.

    Returns:
    - The computed terms for chi.
    """
    # Get the eigenstates and eigenenergies
    eigenenergies, eigenstates = fluxonium.eigensys(evals_count=max_level + 1)

    # Calculate transition frequencies omega_l0 and omega_l1
    omega_10 = eigenenergies[1] - eigenenergies[0]
    omega_rho = ChainMode.get_frequency_from_theory()
    # Initialize the sum for the chi terms
    chi_terms = 0

    # Get the charge operator
    n_operator = fluxonium.n_operator()  # Convert to a full matrix if it is sparse

    # Compute the terms in the square brackets
    for l in range(2, max_level + 1):
        omega_l0 = eigenenergies[l] - eigenenergies[0]
        omega_l1 = eigenenergies[l] - eigenenergies[1]

        # Calculate the matrix elements using the charge operator and eigenstates
        n_0l = eigenstates[:, 0].conj().T @ n_operator @ eigenstates[:, l]
        n_1l = eigenstates[:, 1].conj().T @ n_operator @ eigenstates[:, l]
        n_01 = eigenstates[:, 0].conj().T @ n_operator @ eigenstates[:, 1]

        # Compute the terms and add to the sum
        chi_terms += (np.abs(n_0l)**2) * omega_l0 / (omega_l0**2 - omega_rho**2)
        chi_terms -= (np.abs(n_1l)**2) * omega_l1 / (omega_l1**2 - omega_rho**2)

    # Add the term with the 0-1 transition
    chi_terms += 2 * (np.abs(n_01)**2) * omega_10 / (omega_10**2 - omega_rho**2)
    chi_terms *= ChainMode.coupling_to_fluxonium()**2* 0.5* np.sqrt(ChainMode.josephson_energy_of_chain_junctions/(8*ChainMode.charging_energy_of_mode()))

    return chi_terms

# def get_delta(qubit, evals_count, fc):
#     """
#     Calculate the detuning of transitions from a given frequency.

#     Parameters:
#     qubit (Qubit object): The qubit for which detuning is calculated.
#     evals_count (int): The number of eigenvalues to consider.
#     fc (float): The frequency from which detuning is measured.

#     Returns:
#     ndarray: A 2D array representing the detuning of each transition.
#     """
#     # Get the energy transitions for the qubit
#     transitions = get_transitions(qubit, evals_count)
#     # Subtract the given frequency from each transition to find the detuning
#     return transitions - fc

# def get_condition_dispersif_shift(g_ij, delta):
#     """
#     Calculate the dispersive shift condition for each pair of states.

#     Parameters:
#     g_ij (ndarray): Matrix of coupling strengths between states.
#     delta (ndarray): Matrix of detunings between states.

#     Returns:
#     ndarray: A matrix representing the dispersive shift condition for each pair of states.
#     """
#     # Divide the coupling strengths by the absolute values of the detunings
#     return np.divide(g_ij, np.abs(delta))

# def get_Xllprime(g_ij, delta):
#     """
#     Calculate the cross-Kerr coefficients for each pair of states.

#     Parameters:
#     g_ij (ndarray): Matrix of coupling strengths between states.
#     delta (ndarray): Matrix of detunings between states.

#     Returns:
#     ndarray: A matrix representing the cross-Kerr coefficients for each pair of states.
#     """
#     # Square the coupling strengths and divide by the detunings
#     return np.divide(g_ij**2, delta)

# def get_Xl(Xllprime, index):
#     """
#     Calculate the individual cross-Kerr shift for a specific state.

#     Parameters:
#     Xllprime (ndarray): Matrix of cross-Kerr coefficients between states.
#     index (int): The index of the state for which to calculate the shift.

#     Returns:
#     float: The cross-Kerr shift for the specified state.
#     """
#     # Sum over the rows and columns separately for the specified index
#     X_index_lprime = np.sum(Xllprime[index, :])
#     X_lprime_index = np.sum(Xllprime[:, index])
#     # Subtract the two sums to get the net shift
#     return X_index_lprime - X_lprime_index