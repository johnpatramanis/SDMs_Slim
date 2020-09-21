
###### Files that run
Bash script, runs the Slim file then the python analysis files, if you wanna run the whole pipeline remember to change the K,sigmas in the python file as well


SMD_analysis.py >>
compare_matrices.py >> both are python files that take the outputted file of slim to make the plots

Parameters.txt: The good old paramters file to run the slim script

test_SDM_vector_small_highParea_nrow_40_ncol_60.txt, etc:  The original matrix of SDM values, only used one of them


###### OUTPUTS


The OUT_SED are the output files from slim, one for each run. Contains all information that we output for every individual and every generation, use this if you wanna run your own plots

The 'X_matrix' files are the output of the python file and contain the accumulated corresponding metric (eg competition,number of individuals) over the course of the simulation, for ALL the generations, for each grid.
I collect them and assign them by row, from left to right to wanted size. I remember in R there was an issue when trying to assing them to a matrix whith a certain package.





###### Plot , have sigma and K used in their title


Populations_size_plot:    depicts population over the course of the simulation

Optimum_Distribution_plots:    depicts the enviroment felt by an individual occupying that grid

Occurance_Map: depicts the average number of individuals found in a grid per generation

Distribution_plot: The same as above, not sure why I generated it twice

Normalised_Occurance : the above number but normalised between 0 and 1

Fitness_Distribution: the average fitness of the individuals per grid per generation. (the chance of a person surviving there per generation, env + competition)

Competition_Distribution: the avergage competition felt per individual per grid per generation.

Differences_Plot: the difference between the normilised occurance plot and the origninal SDM values map  (Occurance - Original Value = Difference)



#### Enviromental_Map_plot: The orignal SDM values plotted





##### Correlation Plots, to find where our model performs the worse compared to wanted results (original map values)

Correlation_Between_differences_and_originaMap: Correlation of differences map and the original map values, (some correlation in places where there should be more poeple than what we have in our results)

Correlation_Between_Differences_Competition: Same as above but comparing differences and competition. (No corr)

Correlation_Between_Differences_Fitness: Same as above but comparing differences and fitness. (No corr)
