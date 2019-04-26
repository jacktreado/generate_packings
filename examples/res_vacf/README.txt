running instructions for run_res_velocity_md.sh. The code will take an initially dilute packing, compress it to jamming, and then compress a bit more, energy minimize, and run MD for a set number of steps and ouput velocity information

STEPS:

1. Edit line 9 & 17 to define the path on your machine for the workingdir (directory with generate_packings code) and outdir (directory for program output files)

2. Look at INPUT PARAMETERS section, these will be input to the script at run time and will then be passed into the compiled code. These are:

	N 			: number of residues
	phi0 		: initial packing fraction for packing simulation
	dphi		: compression amount
	dphiDM		: compression amount AFTER jammed packing found
	vacfNT 		: number of time steps to run MD
	vacfT0	 	: temperature of MD run (contact breaking ~ 1e-14)
	tsave	 	: number of steps between velocity saves during MD
	input_str	: path to string with dilute particle data

	NOTE: inputs powers for vacfNT and tsave, the code will then set these variables equal to 2 raised to their input values (for better FFT performance
	)

2. Run script with bash, for example:
	bash run_res_velocity_md.sh 8 0.1 0.001 1e-6 14 1e-20 6 ~/_pv/cluster/rigidbody/io/res_input_N8_seed1.dat 

