Here we share simulation results comparing three Particle-In-Cell (PIC) codes: 
- EPOCH (using the normal random seed and 20 additional runs with different seeds),
- LSP (using both explicit and implicit modes labeled LSP_E and LSP_I), 
- WarpX

Notes:
- Since the previous upload, the WarpX particle energy distributions have been updated to include the virtual (py) momenta
- LSP makes calculations using a virtual dimension of 1 cm while the other codes use 1 meter
- See the corresponding manuscript for full problem description and context

Explanation of Folders/Files:
-> Energy_Diagnostics
	- This directory has files containing field and particle energy information for each code
	- Various file formats are used, see headers for information on data format
	- As mentioned above LSP uses a different convention for the virtual dimension size, so the outputs must be multiplied by 100 for comparison 
	
-> Particle_Spectra
	- This directory contains calculated energy distributions of forward going particles for the various runs. 
	- The file number multiplied by 10 fs gives the simulation time (e.g. 30_Spectrum.dat is 300 fs)
	- In the files, the first column is the energy bin (MeV), the second is total ion charge (C), and third is electron charge (C)
	- As mentioned above LSP uses a different convention for the virtual dimension size, so the outputs must be multiplied by 100 for comparison
	
->Simulation_Input_Files
	- EPOCH and WarpX simulation input decks (input files) are included here 
	- LSP is a commercial code. If you have access to the code and are interested in the input files, please contact us. 
	
-> EPOCH_20_Runs_Density_Map
	- An array of the simulation space denoting whether at least one of the 20 additional EPOCH simulations had non-zero ion density for the corresponding plot step (every 10 fs)
	