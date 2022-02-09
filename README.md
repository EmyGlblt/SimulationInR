# SimulationInR
This repository displays functions for ecological datasets simulations in R. 

## Table of contents
* [Main functions](#Main-functions)
* [Outreach](#Outreach)
* [License](#License)
	
## Main functions
Here are the main function created for simulation and what they can do:
  ### Sample data from a species true intensity:
   * ```simpo``` allows to simulate biased presence-only patterns from an underlying species intensity
   * ```simocc``` allows to simulate biased occupancy data from an underlying species intensity
	
  ### Observer recording and accuracy:
   * ```misidpoints``` simulate misidentification from a supplied point pattern coming from simpo output (simpoRes) 
   * ```multiplepo``` generate multiple species point patterns
    
  ### Data-dynamic:
   * ```deathpoints``` simulate the death of point randomly (random = TRUE) or using habitat suitability (random = TRUE) from a supplied output from simpo, simpoRes 
   * ```birthpoints``` simulate birth (birth=TRUE) and colonization processes (birth=FALSE)
   * ```movepoints``` simulate point movements randomly (mov_type = "random"), depending on habitat suitability (mov_type = "intensity") or a environmental covariate in particular (mov_type = "center").
    
## Outreach
This work is part of my PhD thesis here and has ben presented to the online R in Ecology conference [here](https://github.com/EmyGlblt/SimulationInR/blob/main/Simulation_EcologyinR_EmyGuilbault.pdf)

 
## License
The license for this work is available [here](https://github.com/EmyGlblt/SimulationInR/blob/main/License.txt)
