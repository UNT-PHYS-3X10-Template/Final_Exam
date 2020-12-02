# Final Exam Take-Home Project

In this project you will use a classical molecular dynamics program to study the kinetic energy distribution of a Lennard-Jones gas. 

The project is organized into five tasks, but the core of the program is already implemented.  

The key physical insights that we are going to use to check the quality of our simulation results are basic concepts from statistical thermodynamics, which I summarize in the following:
1. In a statistical system, we expect energy to be equally distributed among equivalent degrees of freedom. The equipartition theorem states that each degree of freedom that contributes quadratically to the total energy of the system should have an average energy of <a href="https://www.codecogs.com/eqnedit.php?latex=\frac{1}{2}k_BT" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{1}{2}k_BT" title="\frac{1}{2}k_BT" /></a>, where <a href="https://www.codecogs.com/eqnedit.php?latex=T" target="_blank"><img src="https://latex.codecogs.com/gif.latex?T" title="T" /></a> is the temperature of the system (in K) and <a href="https://www.codecogs.com/eqnedit.php?latex=k_B" target="_blank"><img src="https://latex.codecogs.com/gif.latex?k_B" title="k_B" /></a> is known as the Boltzmann constant and has a value of 1.38e-23 J/K or 3.165e-6 a.u./K. If you have N particles in D dimensions, the total kinetic energy is thus going to be equal to <a href="https://www.codecogs.com/eqnedit.php?latex=\frac{DN}{2}k_BT" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{DN}{2}k_BT" title="\frac{DN}{2}k_BT" /></a>.
2. The distribution of the speed of the classical distinguishable particles will follow the Maxwell-Boltzmann distribution, which is a close relative of the Gaussian (normal) distribution function. Namely, in three dimensions, the expected distribution function for the speeds should be <a href="https://www.codecogs.com/eqnedit.php?latex=f^{MB}(v)=\left(\frac{m}{2\pi&space;k_BT}&space;\right&space;)^{3/2}4\pi&space;v^2&space;e^{-\frac{mv^2}{k_BT}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f^{MB}(v)=\left(\frac{m}{2\pi&space;k_BT}&space;\right&space;)^{3/2}4\pi&space;v^2&space;e^{-\frac{mv^2}{k_BT}}" title="f^{MB}(v)=\left(\frac{m}{2\pi k_BT} \right )^{3/2}4\pi v^2 e^{-\frac{mv^2}{k_BT}}" /></a>.

TASKS:
1. Initialize the size of the box, the number of particles, and their initial velocities. Note that velocities can be initialized 
   at random, but the velocity of the center of mass of the system needs to be equal to zero.
2. Optimize the timestep of the simuation so that the total energy of the system is conserved
3. Collect the speed and kinetic energy of the particles in time and perform a statistical analysis of the data:
   - What is the mean kinetic energy of each particle?
   - Are there significant differences between different particles? 
   - Plot the frequency distribution of the speed of the particles and compare it to the Maxwell-Boltzmann distribution function.
4. Collect the total kinetic energy of the system and perform a statistical analysis of the data:
   - What is the temperature of the system? Report your estimate with a 99% confidence interval.
   - Plot the frequency distribution of the total kinetic energy of the system. What distribution function do you expect?

OPTIONAL TASKS:
1. Repeat the analysis above for a larger number of particles. How are the results different? 
