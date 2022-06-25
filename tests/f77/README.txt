1. Contents:
The program is in one Fortran source code file, cg.f, that needs to be compiled before usage (see section 2). It is compatible with Fortran77 and newer Fortran standards. Section 3 describes using the program, section 4 input files, section 5 output files and section 6 is the MIT license. This archive contains also 3 different usage examples, each in a different directory:

example1 - Go model simulation of folding a structured globular protein, ubiquitin
example2 - PID model simulation of Atomic Force Microscope (AFM) stretching of an intrinsically disordered protein (IDP), alpha-synuclein
example3 - quasi-adiabatic model simulation of deforming a simulation box filled with partially structured gluten proteins

If you are using this code in a Code Ocean capsule, you can run the examples 1,2 and 3 by running run.sh, run2.sh and run3.sh respectively. run.sh is set as the file to run by default.

===========================================================================


2. Compilation examples with the freely available gfortran compiler:

a) standard compilation:
gfortran cg.f -o cg

b) for debugging purposes:
gfortran -fcheck=all -Og cg.f -o cg

c) for faster execution:
gfortran -O3 cg.f -o cg

d) for the fastest, but less precise calculations:
gfortran -O3 -ffast-math -march=native cg.f -o cg

e) for multi-core calculations with openMP (tested on a Linux machine with 4 cores):
ulimit -s unlimited
gfortran -fopenmp cg.f -o pid
env OMP_NUM_THREADS=4 ./cg inputfile

In the case e) the "ulimit -s unlimited" command increases the stack size so that each thread can make a private copy of a reduced array, which uses a lot of memory (and may cause segmentation fault even for a small protein because array sizes are static). This array reduction makes the multi-core parallelization inefficient for more than 8 cores.

===========================================================================


3. Usage (all the necessary files should be in the same directory):
./cg inputfile

The program takes at most one argument, the name of the inputfile. If no argument is given, the program uses the default values for all variables. The program takes only the inputfile as an argument, names of the other files need to be declared as variables in the inputfile or in the program itself.

Each line of the inputfile should contain the name of the variable and then its value, separated by a single space. Logical variables start with letter l and can be set to T (true) or F (false). Some of them are defined in section 4a.

The program can accept either a protein sequence file (with format defined in section 4c) or a PDB file (with no missing residues). Their names are stored in variables seqfile and pdbfile, respectively. Only one is needed for the simulation. The total number of residues is limited to 10000. If you want to study bigger systems, please replace 10000 with a bigger number in all occurences of "parameter(len=10000)" in the cg.f file.

The program can perform simulations with different kinds of potential:
- structure-based Go-model potential. 
  To use it, set the variable lpdb T (and the contact map will be constructed automatically from the PDB file) or specify the contact map manually in the seqfile
- custom potential for intrinsically disordered proteins.
  To use it, set the variable lcpot T. It may use either a pdbfile or a seqfile. If variables lcpot and lpdb are both set to T, the custom potential will be used alongside the structure-based potential. To use a new version of the custom potential that uses pseudo-improper-dihedral angles, set the variable lpid T (and keep lcpot T).

There are also different potentials for the backbone stiffness of the chain:
- chirality potential (only for structure-based simulations), set by variable lchiral
- bond and dihedral angles potential, set by langle and ldih variables. It has 3 different modes:
a) structure-based (with minima corresponding to the native structure), set by lpdb
b) statistical, based on random-coil library (for IDPs), set by lcoilang
c) user-defined, based on the file tab_ene.txt, set by lenetab
If variable lcpot is set to T, a parameter file (with name stored in variable paramfile) is needed to provide values of the potential parameters. The default paramfile parameters.txt is included in this archive. Its format is defined in section 4d.

The PDB file or the seqfile may contain one or more protein chains (DNA chains in PDB files are ignored). Both single- and multi- chain simulations can be performed using one of many possible boundary conditions:
- infinite space (there are no boundaries and no simulation box, just the protein in implicit solvent). Suitable for single-chain or strongly bound multi-chain complex simulations
- cuboid simulation box, suitable for multi-chain simulations. Different boundary conditions may be used for the walls perpendicular to the Z axis ("in Z direction") and for the walls perpendicular to X and Y axis ("in X and Y directions"). These include flat repulsive wall, flat attractive wall, wall made of face centered cubic (FCC)-packed beads or periodic boundary conditions. Variables controlling them are described in section 4a.

The program produces several files with names stored in the following variables:
outfile - general output of the program, saved every kwrite (in tau time units)
savfile - PDB structures produced by the simulation, saved every ksave (in tau time units)
mapfile - contact map produced by the simulation, saved every ksave (in tau time units)
If you do not want to use the separate name for each of the files, you can specify in the variable "file" a name that will be used for all 3 of them. Then the resulting files will have extensions .out, .pdb and .map. If this option is specified, the number of characters in the "file" variable must be specified in variable klentstr.
The format of all three files is defined in section 5.

===========================================================================


4. Detailed dscription of all the input files that may be used in the simulation:

a) inputfile - list of variables characterizing the simulation. Every line contains one variable name and its value (separated by a single space). Below is the list of more important variables (for the full list and the default values see lines 73 and below in the source code):
	klenstr - length of the output file name (without dot and the extension)
	file - output file name (if klenstr is set to 4 and file is set to lorem, files lore.out and lore.pdb will be generated)
	seqfile - filename of the sequence to be simulated (if there is no input PDB file)
	lparam - logical variable (T or F) that tells if external parameter file is needed
	paramfile - filename of the parameter file (if provided in the inputfile, it automatically sets lparam T)
	lwritemap - logical variable that tells if the output contact map file should be printed during simulation
	mapfile - name of the optional output file with contact map (default is the same as the "file" variable, but with .map extension)
	lwritego - logical, if set to T, no simulations are performed and only a simple Go-model map is printed
	lwritexyz - logical variable that tells if the output conformations should be printed in XYZ format instead of the default PDB format
	iseed - seed for the random number generator (the same seed on the same machine should produce identical simulation)
	ntraj - number of trajectories (independent simulations) to be done
	mstep - time of simulation (in tau time units)
	lcpot - logical, if the custom potential for IDPs (with parameters from paramfile) is turned on
	potcoeff - coefficient multiplying the custom potential (default is 1)
	temp - simulation temperature (in epsilon units, where Boltzmann constant k=1. "Room temperature" varies from 0.3 to 0.7, depending on the model)
	cut - distance for the hard-core repulsion (in Angstroem)
	rcut - potential cutoff (all potentials are set to 0 at this distance. Default is 18 Angstroem)
	verlcut - distance added to Verlet list cutoff (no physical meaning, 10.0 Angstroem is good for performance reasons). Total cutoff is rcut+verlcut
	lpdb - logical, tells if a PDB structure should be used (and Go-model contact map constructed from it)
	lstartpdb - logical, tells if the simulation should start from a given PDB structure (regardless of the model)
	pdbfile - PDB file with the protein to simulate
	lunwrap - logical, allows to load a PDB file saved as one periodic image (with periodic boundary conditions)
	lallatom - logical, tells to construct Go-model based on the Overlap Criterion for heavy atoms (if set to F, the alternative is using Calpha-Calpha distance with cutoff dnaver)
	dnaver - cutoff distance for defining native contacts from an initial pdb file (in Angstroem. Default is 0, which means no native contacts)
	lmedian - logical, for simulations of folding, computes the median time of reaching the native conformation (in tau time units)
	lconftm - logical, for calculating the mean time of formation of each contact (may be used independently of lmedian, time is in tau units)
	lnatend - logical, tells if the trajectory should end when the native conformation is reached
	lunfold - logical, for calculating the mean time of breaking of each contact (may be used independently of lmedian, time is in tau units)
	lvelo - logical, for simulations of pulling protein termini with constant speed of the Atomic Force Microscope tip
	velo - speed of AFM pulling (in angstroem/tau units)
	HH1 - spring constant for pulling AFM springs and attaching beads harmonically to the wall (in epsilon/angstroem^2 units)
	H1 - spring constant for springs connecting neighbour residues in the chain (in epsilon/angstroem^2 units)
	lforce - logical, for simulations of pulling protein termini with constant force by the AFM tip
	screend - electrostatic screening distance
	coul - electrostatic interactions constant (85 is correct for the Tozzini model, for simple D-H electrostatics in water use 2.63, and 210 for vacuum)
	lecperm - logical, T means using constant relative permittivity (simple D-H), F means using Tozzini electrostatics
	lrestart - logical, tells if the simulations should be restarted from the restart file rstfile (default is F)
	rstfile - name of the file used to restart the simulation (the inputfile should be the same as the original, but with two lines setting lrestart T and specifying rstfile name added)
	krst - how often to save the restart file (in tau units; 0 - don't save the restart file)
	ldelrst - if the old restart file should be deleted when the new one is saved
	lcleanrst - rebuild the contact map every krst tau (useful only for the quasi-adiabatic model)
	lmj - if the paramfile should include amplitudes of ss interactions between amino acids (e.g. Miyazawa-Jernigan)
	epsbb - amplitude of backbone-backbone interactions (in epsilon units)
	lsink - logical, for the sink-shaped L-J interaction potentials with a flat region near the minimum
	langle - logical, for the backbone stiffness potential based on bond and/or dihedral angles
	ldih - logical, for dihedral angle potential (works only if langle is T). Named "ldi" in the code
	lcoilang - logical, for using a knowledge-based backbone stiffness angle potential from random coil library (it is read from the paramfile)
	lenetab - use energy values of the bond and dihedral potential from a table called tab_ene.txt
	lchiral - logical, tells if the chirality potential should be used for backbone stiffness (requires pdb structure)
	lmass - logical variable that tells if different amino acid masses should be taken into account (if set to F, each residue has the mean amino acid mass)
	kteql - system equilibration time (in tau time units)
	ksave - how often the output PDB file (and optionally the contact map) should be saved (in tau)
	kksave - ksave value applied before the end of equilibration (if you want to save structures more or less often during equilibration)
	kwrite - how often to save simulation results in the .out file (in tau)
	lmrs - logical, if the morse potential for disulfide bonds is on
	ldynss - logical, controls dynamic disulfide bond formation
	cntfct - determines at what distance the contact is broken (in units of sigma from the L-J potential)
	lpid - logical, for using analytical, pseudo-improper-dihedral potential
	lwall - if there is a simulation box (if F, there are no boundary conditions)
	lfcc - if the walls in the Z direction should be made from 2 FCC layers of beads (otherwise walls are continous)
	lwalls - if there should be repulsive walls in X and Y directions
	lnowal - if there should be no solid (repulsive or attractive) walls in Z direction
	lpbc - logical, for periodic boundary conditions (you need to specify lwall T, lwalls F and lnowal T to get p.b.c. in all 3 directions)
	lcpb - logical, for using periodic boundary conditions even during generation of initial random conformations (the box size does not change even if the random conformation crosses the wall)
	ldens - the system will shrink to reach target density (useful only if lwall T)
	loscillate - for simulations of oscillating deformation of the simulation box (useful only if lwall T)
	lshear - for a shearing oscillations simulation (useful only for lwall T and loscillate T)
	lminfbox - if the box (after shrinking to target density) should be extended in the Z direction until the response force on the Z walls is 0
	lconstvol - keep constant volume during box deformations
	lpullfin - extend the box in the Z direction after oscillations
	densvelo - speed of shrinking (in angstroem/tau units)
	tdens - target density in residues/cubic Angstroem
	sdens - initial density in residues/cubic Angstroem
	wallmindist - minimum of the wall potential
	ktrest - time (in tau) of "resting" (equilibrium simulation) after shrinking (it can be even bigger than mstep, then the system will be "resting" for the rest of the trajectory)
	kconnecttime - when to make the Z walls attractive (1 - just after shrining, 3 - after resting, 5 - after extending the box to the maximal amplitude of oscillations, 7 - after oscillations, 9 - never), only odd values allowed

b) pdbfile - standard PDB file without missing residues (if more than one structure is provided, the first one will be used)
	
c) seqfile - sequence file, consists of the following lines:
- the first line may start with "screend". If it does, after a single space the D-H screening length should be given. The "screend" line is optional, but should be at the beginning of the file).
- the second (or first, if screend is not set in the seqfile) line consists of the number of chains in the seqfile
  Then each chain has the following lines:
  - number of residues in the chain
  - chain sequence in the one-letter format
  - number of contact map files used for the chain
  - names of the contact map files, each in a separate line

The input contact map format (used only if there is no input pdb file provided) is the following:
- the first line is the shift of the reading frame: if it is set to 0, number "1" in the contact map corresponds to the first residue in the chain. If the shift is set to e.g. 5, number "1" corresponds to the 7th residue in the chain.
- the second line is the number of residues considered in the contact map
- the third line is the number of contacts, n
- each of the following n lines contains one contact, each has 3 columns:
-- the first two columns are the numbers of two residues in a contact
-- the 3rd column is the distance between the residues in the units of 5 Angstroem (if the distance is 6 Angstroem, it is written as 1.2 in the contact map)
- then the bond and dihedral angles for each residue need to be provided (in radians), one line for each residue. The terminal angles are undefined, but they also need to be provided (e.g. as zeroes).

d) paramfile - first 22 lines specify the backbone stiffness parameters. 
Lines 23-27 introduce amino acid specificity:
Line 24 - 0 for GLY/PRO, 1 for hydrophobic, 2 for polar, 4 and 5 for charged
Line 25 - coordination number
Line 26 - hydrophobic coordination number
Line 27 - polar coordination number
Amino acid radii are not used, unless variable lradii is set to T
The last part of the paramfile are the potential minima for the pairwise sidechain-sidechain distances for each of 210 possible amino acid pairings. 
If lmj is set to T, a 4th column with the amplitude of the interaction for each pair should be present (in the archive we include parameters.txt without the interaction matrix, parametersMJ96.txt with the Miyazawa-Jernigan matrix from 1996, parametersMDCG.txt and parametersMD01.txt with the MDCG matrix rescaled by 1 or 0.1, respectively.

f) tab_ene.txt - tabularized values of the energy associated with the backbone stiffness bond and dihedral angles (in epsilon units). Each row corresponds to 0.01 degree (around 0.000174532 radians), each column to 9 possible GLY/PRO/X combinations. First 9 columns are for the bond angle potential, last 9 for the dihedral angle potential

===========================================================================


5. Output files (if the "file" variable is set to "name", the filenames will be those in parentheses):

a) outfile (name.out) - simulation output in columns. Description of columns (from left to right) below:
	1. TIME - simulation time (in tau units)
	2. EPOT - potential energy (in epsilon units)
	3. ETOT - total energy (in epsilon units)
	4. ICN - number of native contacts present in a given time
	5-7. B1-B2 B1-S2 S1-S2 - number of backbone-backbone, backbone sidechain and sidechain-sidechain contacts between chains
	8-11. B1-B1 B1-S1 S1-S1 - number of backbone-backbone, backbone sidechain and sidechain-sidechain contacts within one chain
	12. RG - gyration radius
	13. L - end-to-end distance
	14. RMSD - root mean square deviation from the starting structure
	15. NCORD - coordination number
	16. W - asphericity parameter
	17. CORDR - contact order
	18. KNOTS - first residue in the chain that is part of a knot (0 if there is no knot)
	19. KNOTE - last residue in the chain that is part of a knot (0 if there is no knot)
Some of the columns may not be present in some types of simulations (e.g. there will be no end-to-end distance column in a simulation of many chains), and some other columns may be present (for example the response force in simulations that involve deforming the system or the state of the simulation)

b) savfile (name.pdb) - output structure saved in PDB format:
	first line with the CRYST record describes the dimensions of the simulation box (0 for infinite box)
	under every chain there is a REMARK line containing details about that chain (RG, L, W, KNOTS, KNOTE etc.)
	under every chain there is a HETATM line named COG that represents its Center Of Gravity
	If periodic boundary conditions are used, the PDB file contains original, unwrapped coordinates that may extend over many periodic images. In VMD, you need to use "pbc wrap -all" command to put all the residues in one periodic image
	
c) mapfile (name.map) - contact map in format:
	# number of contacts - number of contacts in a given map
	lines showing contacts and backbone stiffness angles in format (columns described from left to right):
	for contacts:
		- K (letter indicating that this line describes a contact)
		- number of a contact (1,2,3...)
		- time of recording the contact map (in tau)
		- number of the 1st residue making that contact
		- number of the 2nd residue making that contact
		- unique id computed for that particular pair of residues
		- type of contact (4 - backbone-backbone, 5 - sidechain-sidechain, 6 or 7 - backbone-sidechain)
	for angles:
		- A (letter indicating this line describes an angle)
		- number of the angle - for a chain with length N there are N-3 dihedral angles and N-2 bond angles (convention: the "first" angle is not printed)
		- time of recording the map (in tau)
		- bond angle with a given number (in radians)
		- dihedral angle with a given number (in radians)

===========================================================================


6. Copyright 2020 Institute of Physics, Polish Academy of Sciences

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
