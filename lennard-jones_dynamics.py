"""
Final Exam Take-Home Project

The following is a program to simulate the motion of objects interacting via a Lennar-Jones potential. The Velocity Verlet algorithm is adopted to 
integrate the equations of motion. 
 
Authors: O. Andreussi and STUDENT
"""
# Import useful modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
# import vpython as vp
# Define basic functions for potential energies and forces
def distance(r1,r2):
    dr=r1-r2
    dr=dr-np.rint(dr/cell)*cell
    return np.sqrt(sum(dr**2)),dr
def apply_pbc(r1):
    r1 = r1 - np.floor(r1/cell)*cell
    return r1
def ulj(x,epsilon=4*4.e-4,sigma=6.):
    """Compute Lennard-Jones potential as a function of inter-atomic distance"""
    return epsilon*((sigma/x)**12-(sigma/x)**6)
def flj(x,epsilon=4*4.e-4,sigma=6.):
    """Compute Lennard-Jones force as a function of inter-atomic distance"""
    return epsilon*6/x*(2*(sigma/x)**12-(sigma/x)**6)
# Define basic functions to compute energy and forces for a collection of objects
def calc_kinetic(masses, velocities):
    """
    Function to compute the total kinetic energy of the system, given the 
    masses and velocities (speeds) of the objects
    
    This function takes in input two numpy arrays and returns a single float number
    """
    speed2 = np.sum(velocities**2, axis=1)
    kinetic_energy=0.5*masses*speed2
    return sum(kinetic_energy)
def calc_forces(positions, f, u):
    """
    Function to compute the total force acting on each object, given the positions
    of the objects and the expression of a central force, i.e. a force that only 
    depends on the pair-wise separation between objects. The function also computes
    the total potential energy of the system.
    
    This function takes in input a numpy array with the x coordinates of n objects
    and a function f that computes a central force as a function of the distance
    between pairs of objects. 
    
    This functions returns in output a numpy arrays with the total force acting 
    on each of the objects.
    """
    nobjects=len(positions)
    forces=np.zeros(positions.shape)
    potential_energy = 0.
    for i in range(nobjects):
        for j in range(i+1,nobjects):
            dist,dr=distance(positions[j,:],positions[i,:])
            potential_energy+=u(dist)
            ftmp=f(dist,epsilon,sigma)*dr/dist
            forces[i,:]-=ftmp[:]
            forces[j,:]+=ftmp[:]
    return potential_energy, forces
# Define functions to initialize the configuration of the system
def init_positions(nobjects,cell):
    """
    Function to initialize the positions of the objects by placing on a regular grid
    
    This function takes in input the number of objects and the cell sizes 
    
    This function returns an array of positions with nobjects * ndim values
    where ndim is the number of elements in the cell array
    """
    ndim = cell.size
    ngrid = int(np.ceil(nobjects**(1./ndim)))
    gridpoints = np.zeros((ngrid**ndim,ndim))
    if ndim == 1 : 
        mesh = np.meshgrid(np.arange(0,cell[0],cell[0]/ngrid))
        gridpoints[:,0]=mesh[0].reshape(-1)
    elif ndim == 2 :
        mesh = np.meshgrid(np.arange(0,cell[0],cell[0]/ngrid),np.arange(0,cell[1],cell[1]/ngrid))
        gridpoints[:,0]=mesh[0].reshape(-1)
        gridpoints[:,1]=mesh[1].reshape(-1)
    elif ndim == 3 :
        mesh = np.meshgrid(np.arange(0,cell[0],cell[0]/ngrid),np.arange(0,cell[1],cell[1]/ngrid),np.arange(0,cell[2],cell[2]/ngrid))
        gridpoints[:,0]=mesh[0].reshape(-1)
        gridpoints[:,1]=mesh[1].reshape(-1)
        gridpoints[:,2]=mesh[2].reshape(-1)
    arr=np.arange(ngrid**ndim)
    np.random.shuffle(arr)
    positions=np.zeros((nobjects,ndim))
    for i in range(nobjects):
        positions[i,:]=gridpoints[arr[i],:]
    return positions   
def init_velocities(nobjects,cell,vmax):
    """
    Function to initialize the positions of the objects by placing on a regular grid
    
    This function takes in input the number of objects and the cell sizes 
    
    This function returns an array of velocities with nobjects * ndim values
    where ndim is the number of elements in the cell array
    
    Be sure to remove the velocity of the center of mass of the system
    """
    ndim = cell.size
    # TASK 1 : initialize the velocities with random numbers from -vmax to vmax
    velocities = 
    # TASK 1 : compute the velocity of the center of mass
    vcm=
    # TASK 1 : shift the velocities so that the center of mass does not move
    
    return velocities 
# Define functions to update kinematic properties of objects using Velocity-Verlet
def update_positions(dt,masses,positions,velocities,forces):
    """
    Function to update the positions of the objects according to Newton's laws
    assuming that the force in the time interval dt is constant. 
    
    This function takes in input the masses, positions, velocities, and forces
    of the objects at time t.
    
    This function modifies the positions to their values at time t+dt.
    """
    nobjects=len(masses)
    for i in range(nobjects):
        positions[i,:]=positions[i,:]+velocities[i,:]*dt+0.5*forces[i,:]/masses[i]*dt**2
        positions[i,:]=apply_pbc(positions[i,:])
    return None
def update_velocities(dt,masses,velocities,forces,newforces):
    """
    Function to update the velocities of the objects using Velocity Verlet algorithm
    
    This function takes in input the masses, velocities, and forces at time t,
    as well as the forces at time t+dt. 
    
    This function modifies the velocities to their values at time t+dt
    """
    nobjects=len(masses)
    for i in range(nobjects):
        velocities[i,:]=velocities[i,:]+0.5*(forces[i,:]+newforces[i,:])/masses[i]*dt
    return None
# setup and update vpython visualization
def init_visualization(cell,positions):
    """
    Function to setup the VPython visualization of the system
    """
    vp.scene.caption = """
Right button drag or Ctrl-drag to rotate "camera" to view scene.
To zoom, drag with middle button or Alt/Option depressed, or use scroll wheel.
On a two-button mouse, middle is left + right.
Shift-drag to pan left/right and up/down.
Touch screen: pinch/extend to zoom, swipe or two-finger rotate."""
    ndim = cell.size
    local_cell = np.zeros(3)
    thk = 0.3
    side = 15.0
    s2 = 2*side - thk
    s3 = 2*side + thk
    if ndim == 1 :
        local_cell[0] = cell[0]
        local_cell[1] = cell[0]/4
        local_cell[2] = cell[0]/4
        wallR = vp.box(pos=vp.vector(local_cell[0], 0, 0), size=vp.vector(thk, local_cell[1], local_cell[2]),  color = vp.color.red)
        wallL = vp.box(pos=vp.vector(0, 0, 0), size=vp.vector(thk, local_cell[1], local_cell[2]),  color = vp.color.red)    
    elif ndim == 2 :
        local_cell[0] = cell[0]
        local_cell[1] = cell[1]
        local_cell[2] = (cell[0]+cell[1])/8
        wallR = vp.box(pos=vp.vector(local_cell[0], local_cell[1]/2, 0), size=vp.vector(thk, local_cell[1], local_cell[2]),  color = vp.color.red)
        wallL = vp.box(pos=vp.vector(0, local_cell[1]/2, 0), size=vp.vector(thk, local_cell[1], local_cell[2]),  color = vp.color.red)    
        wallB = vp.box(pos=vp.vector(local_cell[0]/2, 0, 0), size=vp.vector(local_cell[0], thk, local_cell[2]),  color = vp.color.blue)
        wallT = vp.box(pos=vp.vector(local_cell[0]/2, local_cell[1], 0), size=vp.vector(local_cell[0], thk, local_cell[2]),  color = vp.color.blue)
    elif ndim == 3 :
        local_cell[:] = cell[:]
        wallR = vp.box(pos=vp.vector(local_cell[0], local_cell[1]/2, local_cell[2]/2), size=vp.vector(thk, local_cell[1], local_cell[2]),  color = vp.color.red)
        wallL = vp.box(pos=vp.vector(0, local_cell[1]/2, local_cell[2]/2), size=vp.vector(thk, local_cell[1], local_cell[2]),  color = vp.color.red)    
        wallB = vp.box(pos=vp.vector(local_cell[0]/2, 0, local_cell[2]/2), size=vp.vector(local_cell[0], thk, local_cell[2]),  color = vp.color.blue)
        wallT = vp.box(pos=vp.vector(local_cell[0]/2, local_cell[1], local_cell[2]/2), size=vp.vector(local_cell[0], thk, local_cell[2]),  color = vp.color.blue)
        wallBK = vp.box(pos=vp.vector(local_cell[0]/2, local_cell[1]/2, 0), size=vp.vector(local_cell[0], local_cell[1], thk), color = vp.color.gray(0.7))
    vp.scene.center = 0.5*vp.vector(local_cell[0],local_cell[1],local_cell[2])
    spheres=[]
    for i in range(nobjects):
        if ndim == 1 :
            local_position=vp.vector(positions[i,0],0,0)
        elif ndim == 2 :
            local_position=vp.vector(positions[i,0],positions[i,1],0)
        elif ndim == 3 :
            local_position=vp.vector(positions[i,0],positions[i,1],positions[i,2])
        spheres.append(vp.sphere(pos=local_position))
    return spheres
def update_visualization(cell,spheres,positions):
    """
    Function to update the positions of the spheres in the visualization
    """
    ndim = cell.size
    for i in range(len(spheres)):
        if ndim == 1 :
            local_position=vp.vector(positions[i,0],0,0)
        elif ndim == 2 :
            local_position=vp.vector(positions[i,0],positions[i,1],0)
        elif ndim == 3 :
            local_position=vp.vector(positions[i,0],positions[i,1],positions[i,2])
        spheres[i].pos=local_position
    return
#
# Main Program
#
# Define Lennard-Jones parameters for Argon atoms
#
epsilon=4*4.e-4 # (in atomic units) NOTE: added a factor of four with respect to assigments
sigma=6.0 # (in atomic units)
mass=73000. # (in atomic units)
# Define basic setup of the simulation: number of objects, masses, initial positions
# initial velocities, timestep.
#
cell=np.array([20.0, 20., 20.]) # TASK 1: initialize the box (should probably be a 3D box, but you can play with 1D and 2D)
nobjects=8 # TASK 1: select the number of particles
masses=mass*np.ones(nobjects) # Masses of the objects (in atomic units)
positions=init_positions(nobjects,cell) # Positions of the objects (in atomic units)
velocities=init_velocities(nobjects,cell,0.005)# Initial velocities of the two objects (in atomic units)
#
# Setup VPython visualization
#
# spheres=init_visualization(cell,positions)
# 
# Compute the force at the initial step
#
utot,forces=calc_forces(positions,flj,ulj)
#
# During the simulation we will keep track of the following quantities
#
times =[0.]
potential=[utot] # The total potential energy of the system
kinetic=[calc_kinetic(masses,velocities)] # The total kinetic energy of the system
#
# Start the loop for the time integration
#
time=0.
dt=50 # Time step used for the numerical integration (in atomic units) TASK 2: This value needs to be changed
total_time=40000 # Total time simulated (in atomic units)
stride=100 # Save information at every stride iterations
while time<total_time :
#    vp.rate(1000)
    time+=dt # Take one step in time
    update_positions(dt,masses,positions,velocities,forces) # Update positions
#    update_visualization(cell,spheres,positions)
    utot,newforces=calc_forces(positions,flj,ulj) # Compute forces at new positions
    update_velocities(dt,masses,velocities,forces,newforces) # Update velocities
    forces=newforces # Overwrite forces for next step
    #
    # Store the current values of the tracked quantities
    #
    if int(time/dt)%stride == 0 : 
        times.append(time)
        potential.append(utot)
        kinetic.append(calc_kinetic(masses,velocities))
#
# Alternative plots of the two components of the mechanical energy of the system
#
total=[potential[i]+kinetic[i] for i in range(len(potential))]
#
# Plot the energies of the system to check the stability of the dynamics (TASK 2)
#

#
# TASK 3: Analyze the individual kinetic energy of the particles
#

#
# TASK 4: Analyze the total kinetic energy of the system
#
