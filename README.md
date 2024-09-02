# Evolutionary_Rescue_Two_Demes
This repository contains five .cpp files:

(1) The code fixation-min_distance.cpp simulates the fixation process of a single mutation for the case in which the phenotypic optima of demes 1 and 2 are aligned. 

(2) The code fixation-general.cpp simulates the fixation process of a single mutation for the general case, i.e.,  the phenotypic optima of demes 1 and 2 have random orientations. 

(3) The code metapopulation-onestep_min.cpp simulates the rescue process for the case in which the phenotypic optima of demes 1 and 2 are aligned. Only one-step mutations are allowed. 

(4) The code metapopulation-min_distance.cpp simulates the rescue process for the case in which the phenotypic optima of demes 1 and 2 are aligned. Double and multiple mutations are allowed.

(5) The code metapopulation-phenotypic.cpp simulates the rescue process for the general case, i.e.,  the phenotypic optima of demes 1 and 2 have random orientations. Double and multiple mutations are allowed.


All codes above are C++ source code "name.cpp"

Requires the GSL scientific library:
	
Ubuntu systems installation through terminal:
		
    (C++)	~ sudo apt-get install build-essential
		(GSL)	~ sudo apt-get install gsl-bin
		(GSL)	~ sudo apt-get install libgsl-dev

Ubuntu systems compilation through terminal:
	
		~ c++ -O3 name.cpp -o [executable_name] -lm -lgsl -lgslcblas

		After compilation, we will get an executable.

    To run the code, we must provide input data.
                    
    The input data for codes (1) and (2) are:
                      
         i - number of subpopulations (which is equal to two in our simulations)     
         ii- carrying capacity K
         iii - number of traits n
	       iv - migration probability  
	       v - number of independent runs   
	       vi - the strength of selection, which we set at 1
         vii - mean phenotypic effect of the mutant
         viii - stress level, delta_1 (deme1)
         ix - stress level, delta_2 (deme2)
          x - Maximum fitness W_max
        
         Here is an example of how to run the code in an Ubuntu terminal

         ./[executable_name] 2 10000 5 0.001 10000 1. 0.25 0.1 0.2 1.5

        The output of the code is created in a.DAT file. The name of the file includes the values of the parameter 
        used. The first column of the DAT file is delta_1, whereas the second column provides the 
        probability a mutation fixates in demes 1 and 2


     The input data for code (3) are:
                      
         i - number of subpopulations (which is equal to two in our simulations)     
         ii- carrying capacity K
         iii - number of traits n
         iv - mutation probability mu
	       v - migration probability  
	       vi - number of independent runs   
	       vii - the strength of selection, which we set at 1
         viii - mean phenotypic effect of the mutant
         ix - stress level, delta_1 (deme1)
         x - stress level, delta_2 (deme2)
         xi - Maximum fitness W_max
        
         Here is an example of how to run the code in an Ubuntu terminal
  
        ./[executable_name] 2 10000 5 0.005  0.001 10000 1. 0.25 0.1 0.2 1.5

        The output of the code is created in a.DAT file. The name of the file includes the values of the parameter 
        used. The first column of the DAT file is delta_1, the second column provides the 
        probability of extinction, the third column is the probability both demes are rescue

       
    The input data for codes (4) and (5) are:

         i - number of subpopulations (which is equal to two in our simulations)     
         ii- carrying capacity K
         iii - number of traits n
         iv - mutation probability mu
	       v - migration probability  
	       vi - number of independent runs   
	       vii - the strength of selection, which we set at 1
         viii - mean phenotypic effect of the mutant
         ix - stress level, delta_1 (deme1)
         x - stress level, delta_2 (deme2)
         xi - Maximum fitness W_max
         xii - maximum time of simulation (which we assume is 2,000)

          Here is an example of how to run the code in an Ubuntu terminal

          ./[executable_name] 2 10000 5 0.005  0.001 10000 1. 0.25 0.1 0.2 1.5 2000

           The output of the code is created in a.DAT file. The name of the file includes the values of the 
          parameter used. The first column of the DAT file is delta_1, the second column provides the 
        probability of extinction, the third column is the probability both demes are rescue
