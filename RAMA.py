#Import various libraries that are needed for this script to work (generally, we only need to import numpy, Pylians3 
# integration library and a couple of standart python suites, such as sys, os etc.)
import numpy as np
import integration_library as IL

import sys, platform, os
import numpy as np

#Show script logo
file_logo = open('logo.txt', 'r')
file_logo = file_logo.read()
print(file_logo+"\n")

#Define output folders for both neutrino and regular cases respectively
output_nu   = str(input("Directory for reps transfer functions with massive neutrino: "))
output_nonu = str(input("Directory for reps transfer functions without neutrino: "))

#Define folders for matter power spectrum output of reps
output_spectrum_nonu =  str(input("Directory for reps power spectra with massive neutrino: "))
output_spectrum_nu   =  str(input("Directory for reps power spectra with massive neutrino: "))

#Define initial redshift and number of outputs (note that redshifts at which outputs are present should be created using
# numpy.linspace() suite, for example, if outputs are present up to z=49 and there are 500 outputs, use 
# np.linspace(0,49,500) in reps)
zin = int("Initial redshift (must be an integer): ")
output_number = int("Number of reps outputs")
arr = np.linspace(0,zin,output_number)

#Define output file name
output_filename = str("File name for output .txt")

#Import present time values of linear growth rates for both neutrino and regular cosmologies
kh = np.loadtxt(str(output_nu)+"/_rescaled_transfer_z0.0000.txt")[:,0]
Delta_nonu_targetz0 = np.loadtxt(str(output_nu)+"/_rescaled_transfer_z0.0000.txt")[:,7]
Delta_tot_targetz0 = np.loadtxt(str(output_nu)+"/_rescaled_transfer_z0.0000.txt")[:,6]
Delta_nonu_originalz0 = np.loadtxt(str(output_nonu)+"/_rescaled_transfer_z0.0000.txt")[:,7]
Delta_tot_originalz0 = np.loadtxt(str(output_nonu)+"/_rescaled_transfer_z0.0000.txt")[:,6]

#Suite to calculate linear mass fluctuations in real space using Pylians3 integration library
def sigma(k,Pk,R):
    yinit = np.array([0.0], dtype=np.float64)
    eps   = 1e-13  #change this for higher/lower accuracy
    h1    = 1e-12
    hmin  = 0.0
    
    W   = 3.0*(np.sin(k*R) - k*R*np.cos(k*R))/(k*R)**3
    Pk1 = Pk*W**2*k**2/(2.0*np.pi**2)
    
    return np.sqrt(IL.odeint(yinit, 10**(-4),10**(0), eps,
                             h1, hmin, np.log10(k), Pk1,
                             'sigma', verbose=False)[0])

#Import present day Pcb(k) and Pnu(k). It is possible to provide your own CAMB formatted outputs for matter power spectrum 
# if one assumes non-standard cosmology (for example, using MGCAMB or MGCLASS II)
k0, Pk0_nonu = np.loadtxt(str(output_spectrum_nonu)+"/_Pm_rescaled_z0.0000.txt", unpack=True)
k0, Pk0_nu = np.loadtxt(str(output_spectrum_nu)+"/_Pm_rescaled_z0.0000.txt", unpack=True)

#Define target cosmology sigma(R) as a subroutine for the sake of simplicity
def target_sigma(R):
    return sigma(k0,Pk0_nu,R)
        

#Iterate through reps outputs
def derive_sigma(k0):
    
    scaled_redshift = []
    scaled_length = []
    mean_arr = []

    for counter in range(output_number):
        kh = np.loadtxt("output_nu/_rescaled_transfer_z" + str(format(arr[counter], '.4f')) + ".txt")[:,0]
        Delta_nonu_target = np.loadtxt("output_nu/_rescaled_transfer_z" + str(format(arr[counter], '.4f')) + ".txt")[:,7]
        Delta_tot_target = np.loadtxt("output_nu/_rescaled_transfer_z" + str(format(arr[counter], '.4f')) + ".txt")[:,6]
        Delta_nonu_original = np.loadtxt("output_nonu/_rescaled_transfer_z" + str(format(arr[counter], '.4f')) + ".txt")[:,7]
        Delta_tot_original = np.loadtxt("output_nonu/_rescaled_transfer_z" + str(format(arr[counter], '.4f')) + ".txt")[:,6]
    
        pk_target = ((Delta_nonu_target/Delta_nonu_targetz0)**2*Pk0_nu)
        pk_target_tot = ((Delta_tot_target/Delta_tot_targetz0)**2*Pk0_nu)
        pk_original = ((Delta_nonu_original/Delta_nonu_originalz0)**2*Pk0_nonu)
        pk_original_tot = ((Delta_tot_original/Delta_tot_originalz0)**2*Pk0_nonu)
    
        #Array for real-space radial coordinate, that runs from R1=1 to R2=10
        RR = np.linspace(1,10,10)
    
        #Array for scaling, we assume that s in [0.5,1.5]
        arr2=np.linspace(0.5,1.5,100)
        for y in range(len(arr2)):
            mean = 0
            for x in range(10):
                mean+=(np.abs((sigma(k0,pk_original,RR[x]/arr2[y])-target_sigma(RR[x]))/target_sigma(RR[x])*100))
            
            #Calculate mean value for sigma(R) percentage difference over all R
            mean = np.sum(mean)/len(RR)
        
            #Write scaled length, redshift and mean percentage difference to an array
            scaled_redshift.append(arr[counter])
            scaled_length.append(arr2[y])
            mean_arr.append(mean)

    return scaled_redshift, scaled_length, mean_arr

#Create file to write all outputs
def write_output(scaled_redshift,scaled_length,mean_arr,output_filename):
    with open(str(output_filename)+".txt", "a") as outputfile:
        outputfile.write("#zo, s, difference \n")
        for i in range(len(mean_arr)):
            outputfile.write(str(scaled_redshift[i])+","+str(scaled_length[i])+","+str(mean_arr[i])+"\n")
        outputfile.close()

def main():
    scaled_redshift, scaled_length, mean_arr = derive_sigma(k0)
    write_output(scaled_redshift,scaled_length,mean_arr,output_filename)
    
main()
