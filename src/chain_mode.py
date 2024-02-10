"""
chain_mode.py: This module contains the ChainMode class which models the modes in the fluxonium inductive chain, 
based on the paper: 

    - Giovanni Viola and Gianluigi Catelani. "Collective modes in the fluxonium qubit." 
        Physical Review B, vol. 92, no. 224511, 2015. 
        DOI: 10.1103/PhysRevB.92.224511.
"""

import scqubits as scq
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style as mplstyle
import os 

# Get the directory where the current script resides
dir_path = os.path.dirname(os.path.realpath(__file__))
# Construct the full path to the style file
style_path = os.path.join(dir_path, 'mystyle.mplstyle')

mplstyle.use(style_path)
# Define my custom color cycle
custom_colors = ["#264653", "#2a9d8f", "#8ab17d", "#e9c46a", "#f4a261", "#e76f51", "#b43718", "#264653", "#2a9d8f"]
# Set the color cycle
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=custom_colors)


class ChainMode(scq.Oscillator):
    """
    ChainMode is an extension of scq.Oscillator to model a specific quantum chain mode.
    It includes additional parameters and methods relevant to this specific mode.
    """
    # Define allowed additional keys for this class
    ALLOWED_KEYS = {
        'number_of_chain_junctions', 'josephson_energy_of_chain_junctions', 
        'capacitance_to_ground_of_chain_junctions', 'area_of_chain_junctions', 
        'charging_energy_of_fluxonium'
    }

    def __init__(self, chain_mode_frequency, truncated_dim, chain_mode_index, **kwargs):
        """
        Initialize the ChainMode object.

        Parameters:
        - chain_mode_frequency: Frequency of the chain mode. if you want, i prefere to calculate the frequency from the theory.
        - truncated_dim: Dimension used for matrix diagonalization.
        - chain_mode_index: Index of the chain mode.
        - **kwargs: Additional keyword arguments for extended properties. 
				it is recommended to give the values of the following entities :
					- number_of_chain_junctions
					- josephson_energy_of_chain_junctions in GHz
					- capacitance_to_ground_of_chain_junctions in fF
					- area_of_chain_junctions in µm²
					- charging_energy_of_fluxonium in GHz
        """
        # Initialize the superclass with necessary parameters
        super().__init__(E_osc=chain_mode_frequency, truncated_dim=truncated_dim)
        self.chain_mode_index = chain_mode_index

        # Set additional parameters from kwargs
        for key, value in kwargs.items():
            if key in self.ALLOWED_KEYS:
                setattr(self, key, value)
            else:
                setattr(self, key, None)
                print(f'The {key} is set to None.')

        # Update oscillator parameters based on new properties
        self.update_losc()

    def charging_energy_of_chain_junction(self):
        """Calculate the charging energy of the chain junction."""
        return 19.37 / (45 * self.area_of_chain_junctions)

    def inductive_energy_of_chain(self):
        """Calculate the inductive energy of the chain."""
        return self.josephson_energy_of_chain_junctions / self.number_of_chain_junctions

    def charging_energy_of_junction_to_ground(self):
        """Calculate the charging energy of the junction to the ground."""
        return 19.37 / self.capacitance_to_ground_of_chain_junctions

    def charging_energy_of_mode(self):
        """
        Calculate the effective charging energy of the mode.
        This is derived from the charging energy of the chain junction and the junction to ground.
        """
        return (1 / self.charging_energy_of_chain_junction() + 1 / (4 * self.charging_energy_of_junction_to_ground() ) *( 1 / self.sine_mode()**2  + 
        	self.oddity_mode()* ( - 2* self.lambda_parameter() * self.cosine_mode()**2 / (self.number_of_chain_junctions *(self.number_of_chain_junctions-1) * self.sine_mode()**4 )  + 
        		 (self.lambda_parameter() * 2 / (self.number_of_chain_junctions *(self.number_of_chain_junctions-1)) ) **2 * self.summation_for_odd_mode_frequency())))**-1

    def sine_mode(self):
        """Calculate the sine mode of the chain mode index."""
        return np.sin(0.5 * np.pi * self.chain_mode_index / self.number_of_chain_junctions)

    def cosine_mode(self):
        """Calculate the cosine mode of the chain mode index."""
        return np.cos(0.5 * np.pi * self.chain_mode_index / self.number_of_chain_junctions)

    def oddity_mode(self):
        """Calculate the oddity of the chain mode index."""
        return (1 - (-1)**self.chain_mode_index ) / 2  
    def eventy_mode(self):
        """Calculate the oddity of the chain mode index."""
        return (1 + (-1)**self.chain_mode_index ) / 2  

    def lambda_parameter(self):
        """
        Calculate the lambda parameter, which is a ratio involving the total capacitance to ground
        and the capacitance of the fluxonium.
        """
        Ct = 19.37 / self.charging_energy_of_fluxonium + (self.number_of_chain_junctions - 1) * self.capacitance_to_ground_of_chain_junctions
        return (self.number_of_chain_junctions - 1) * self.capacitance_to_ground_of_chain_junctions / Ct

    def summation_for_odd_mode_frequency(self):
        """
        Perform a summation calculation used in the charging energy of mode for odd chain mode frequencies.
        It involves a summation over specific sine and cosine terms of the chain mode indices, excluding the current mode index.
        """
        odd_indices = np.arange(1, self.number_of_chain_junctions, 2)
        odd_indices = odd_indices[odd_indices != self.chain_mode_index]

        sine_terms = np.sin(0.5 * np.pi * odd_indices / self.number_of_chain_junctions)
        cosine_terms = np.cos(0.5 * np.pi * odd_indices / self.number_of_chain_junctions)

        numerator = self.cosine_mode()**2 * cosine_terms**2
        denominator = self.sine_mode()**2 * sine_terms**4 - self.sine_mode()**4 * sine_terms**2

        return np.sum(numerator / denominator)

    def phi_zpf(self):
        """Calculate the zero-point fluctuation of the phase."""
        return (2 * self.charging_energy_of_mode() / self.josephson_energy_of_chain_junctions)**0.25

    def n_zpf(self):
        """Calculate the zero-point fluctuation of the number operator."""
        return (self.josephson_energy_of_chain_junctions / (32 * self.charging_energy_of_mode()))**0.25

    def chain_mode_impedance(self):
        """
        Calculate the impedance of the chain mode.
        Placeholder method - to be implemented if needed.
        """
        print('Method not implemented yet! -_-')
        pass 

    def update_losc(self):
        """
        Update the l_osc parameter based on the current state.
        This should be called whenever any parameters affecting l_osc are changed.
        """
        n_zpf = self.n_zpf()
        self.l_osc = 1 / (np.sqrt(2) * n_zpf)

    def coupling_to_fluxonium(self):
        """
        Calculate the coupling energy to a fluxonium qubit.
        This is based on the charging energy of the fluxonium, mode, and junction to ground,
        as well as the sine and cosine modes.
        """
        return self.eventy_mode()*(4 * self.charging_energy_of_fluxonium * self.charging_energy_of_mode() * 
                self.cosine_mode() / 
                (np.sqrt(2 * self.number_of_chain_junctions) * 
                 self.charging_energy_of_junction_to_ground() * self.sine_mode()**2))

    def get_frequency_from_theory(self):
        """
        Calculate the frequency of the chain mode based on theoretical considerations.
        This involves the square root of the product of 8 times the charging energy of the mode and the Josephson energy.
        """
        return np.sqrt(8 * self.charging_energy_of_mode() * self.josephson_energy_of_chain_junctions)

    def set_frequency_to_theory(self):
        """
        Set the oscillator frequency (E_osc) to the theoretical frequency calculated for the chain mode.
        """
        self.E_osc = self.get_frequency_from_theory()


    def plot_dispersion_and_coupling(self):
        """
        Plot the dispersion relation and coupling strength of the chain modes to the fluxonium.
        Automatically uses the number of chain junctions as the range for the mode index (rho).
        """
        # Generate an array of mode indices
        original_chain_mode_index = self.chain_mode_index
        mode_indices = np.arange(1, self.number_of_chain_junctions + 1)

        # Preallocate arrays for coupling strengths and frequencies
        coupling_strengths = np.zeros_like(mode_indices, dtype=float)
        mode_frequencies = np.zeros_like(mode_indices, dtype=float)

        # Calculate coupling strength and frequency for each mode
        for i, mode_index in enumerate(mode_indices):
            self.chain_mode_index = mode_index
            coupling_strengths[i] = self.coupling_to_fluxonium()
            mode_frequencies[i] = self.get_frequency_from_theory()
        self.chain_mode_index = original_chain_mode_index

        # Plotting setup
        fig, (ax_coupling, ax_frequency) = plt.subplots(1, 2, figsize=(10, 5))

        # Plotting the coupling strength
        ax_coupling.plot(mode_indices, coupling_strengths, 'x')
        ax_coupling.set_yscale('log')
        ax_coupling.set_ylabel(r'$g_\text{coupling}$ (GHz)')
        ax_coupling.set_xlabel(r'Mode Index ($\rho$)')
        ax_coupling.grid(True)

        # Plotting the mode frequency
        ax_frequency.plot(mode_indices, mode_frequencies, 'x')
        ax_frequency.set_ylabel(r'Mode Frequency ($\omega_\rho$) in GHz')
        ax_frequency.set_xlabel(r'Mode Index ($\rho$)')
        ax_frequency.grid(True)

        # Adjust layout
        ax_coupling.set_ylim(1e-8, 1)
        fig.tight_layout()
        plt.show()
        return fig, mode_indices, coupling_strengths, mode_frequencies