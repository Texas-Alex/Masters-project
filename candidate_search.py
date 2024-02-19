#!/usr/bin/env python3
import numpy as np
import sys
import math
import matplotlib.pyplot as plt
import csv
import pyne2001
from astropy import units as u
from astropy.coordinates import SkyCoord
from pygdsm import GlobalSkyModel2016
import yaml

G  = 4.3009172706e-3 	#pc M_0^-1 kms^-2


#Read the data from a .yaml file  ####ADD IN README WHAT THIS LOOKS LIKE#####
def input_vals(filename):
    with open(filename ,'r') as file:
        inputs = yaml.safe_load(file) #dictionary
    return inputs


#Read a gaia_source table and extract the relevant columns as numpy.arrays
#Relevant columns in this case are: Source ID, G-band mean magnitude,
#Parallax, Excess noise and its significance, and the RA/DEC.
def read_gaia(filename):
    try:
        with open(filename, "r") as file: #get gaia headers
            headers = file.readline().strip().split(',')
    except FileNotFoundError:
        print("filename not valid (file not found)")
    
    for i in range(len(headers)): #get the ids of each column
        if headers[i] == 'source_id':
            id_id = i
        if headers[i] == 'phot_g_mean_mag':
            flux_id = i
        if headers[i] == 'parallax':
            parallax_id = i
        if headers[i] == 'astrometric_excess_noise':
            e_id = i
        if headers[i] == 'astrometric_excess_noise_sig':
            sig_id = i
        if headers[i] == 'duplicated_source':
            dup_id = i
        if headers[i] == 'ra':
            ra_id = i
        if headers[i] == 'dec':
            dec_id = i

    

    id, flux,parallax, e, e_sig, dups, ra, dec = np.genfromtxt(filename, unpack = True, usecols=(id_id, flux_id, parallax_id, e_id,sig_id, dup_id, ra_id, dec_id ),delimiter = ',',dtype = str)
    
    return id[1:], flux[1:].astype(float), parallax[1:].astype(float), e[1:].astype(float),e_sig[1:].astype(float), dups[1:], ra[1:], dec[1:]

#Determines the Mass of a star using the Mass-Luminosity
#relation.
def companion_mass(flux, dist):
    L = 4*(math.pi)*(dist**2)*flux
    M = L**(1/3.5)
    return M
    

###The below calculations all make use of the pulsar mass
###function to determine a number of parameters.
###https://www.aanda.org/articles/aa/full_html/2018/04/aa31928-17/aa31928-17.html


#Calculate the projected semi-major orbital axis (ap) for a given
#Gaia object,assuming i=90 degrees and that the "wobble" is not
#zero.
def ap_calculation(mc,mpsr, e, distance, i=90):
    i = math.radians(i)
    ap_sin = (e*distance) * (mc/mpsr)
    return (ap_sin*math.sin(i))


#Calculate the orbital period given the mass of a pulsar and companion,
#and the semi-major axis.
def orbital_period(ap, mc, mpsr):
    top = ap**3 * 4 * (math.pi)**2
    bottom = G * (mc + mpsr)
    return math.sqrt(top/bottom)

#Calculate the mass function based on the companion mass
def mass_function_m(mc,mpsr,i=90):
    i = math.radians(i)
    top = mc * math.sin(i)
    bottom = (mc + mpsr)**2
    return((top)**3/bottom)

#Calculate mass function given the orbital period
def mass_function_pb(ap, pb, i=90):
    i = math.radians(i)
    top = 4 * (math.pi)**2 * (ap * math.sin(i))**3
    bottom = G * pb**2

    return (top/bottom)



def main():

    #Variable definition
    candidates = []
    gaiaSource = sys.argv[1] #filename of the gaia_source file
    inputValues = sys.argv[2] #.yaml file with the required parameters
    gsm = GlobalSkyModel2016()

    id, flux, parallax, e,e_sig, dups, ra, dec  = read_gaia(gaiaSource)
    inputs = input_vals(inputValues)


    #Merge the individual column arrays into a single array.
    merged_gaia = np.column_stack((id, flux ,parallax, e,e_sig, dups, ra, dec ))

    for i in range(len(merged_gaia)):
        source_id = merged_gaia[i][0].astype(float)
        flux_val = merged_gaia[i][1].astype(float)
        parralax_val = merged_gaia[i][2].astype(float)
        e_val = merged_gaia[i][3].astype(float)
        e_sig = merged_gaia[i][4].astype(float)
        duplicated = merged_gaia[i][5]
        RA = merged_gaia[i][6].astype(float)
        DEC = merged_gaia[i][7].astype(float)



        if(i%1000 == 0):
            print("iteration ", i,  "out of ", len(merged_gaia))

        # #Candidates too close to the galactic centre are hard to detect and can be
        # #ignored.
        # if((RA >263 and RA <269)and (DEC >-32 and DEC <-25)):
        #     continue



        dist = 1/(parralax_val)

        #Check if the distance is a positive value and that it is within the galaxy.
        if (dist < 30 and dist > 0):
            
            #Define variables.
            M_c = companion_mass(flux_val,dist)
            ap = ap_calculation(M_c, inputs['M_psr'], e_val,dist)
            pb = orbital_period(ap,M_c, inputs['M_psr'])
            mf = mass_function_m(M_c,inputs['M_psr'])
            mass_funct_diff = 1-(mass_function_m(M_c, inputs['M_psr'])/mass_function_pb(ap,pb))


            #Search based on the values provided in the .yaml file.

            #Not a duplicate source
            if(duplicated != 'true'):
                #Companion mass within reasonable range
                if(M_c <= 200):  

                    #Excess noise parameter check
                    if(e_val >=inputs['ecc'][0] and  e_val<= inputs['ecc'][1]):

                        #Comparing the mass function calculations
                        if(mass_funct_diff <= inputs['percent_err']):

                            #Mass function check
                            if(mf >= inputs['Mass_frac'][0] and mf <= inputs['Mass_frac'][1]):

                                #Orbital period check
                                if(pb >= inputs['Pb'][0] and pb <= inputs['Pb'][1]):

                                    #Determine the Dispersion measure (DM) using the pyne2001 package.
                                    freqs = np.array([150])
                                    c = SkyCoord(dec= DEC*u.degree, ra = RA*u.degree,frame = 'icrs', unit = 'deg')
                                    DM =  pyne2001.get_dm(c.galactic.l.degree, c.galactic.b.degree, dist) #dispersion measure determination

                                    #DM check
                                    if(DM >= inputs['DM'][0] and DM <= inputs['DM'][1]):
                                        T_sky = gsm.get_sky_temperature(c,freqs) #Sky temperature determination using GSM package

                                        #Sky temeperature check
                                        if(T_sky <= inputs["temp"]):
                                            candidates.append([source_id.astype(str),RA, DEC,flux_val,parralax_val,e_val,e_sig, M_c, pb,ap, DM, mf, T_sky])
                                                            #id, ra, dec, flux, parallax, e, companion mass, orbital period, DM, mass function(using masses), sky temp
    
    #Write final candidates to a .csv file
    with open("potential_cands.csv", "w") as f:
        header = ['Gaia source ID', 'RA(deg)', 'DEC(deg)', 'photometric G-band flux (Flux[e^- s^-1])', 'Parallax(mas)', 'astrometric wobble (mas)','wobble signifigance', 'Companion mass (Msun)', 'Orbital period(days)','orbital separation [ap]', 'DM(pc cm^-3)', 'Mass function', 'Sky Temperature(K)']
        writer =csv.writer(f)
        writer.writerow(header)
        for i in range(len(candidates)):
            writer.writerow(candidates[i])
    print("File Written")

    return

if __name__ == "__main__":
    main()