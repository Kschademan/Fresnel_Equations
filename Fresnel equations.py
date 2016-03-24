################################################################################
#Programmer: Kyle Schademan
#Date: March 18, 2016
#Description: This program calculates the S and P polarization trasmitance 
#    through a specificied lens.  This module uses the Fresnel Equations
#    to calculate the polarization
################################################################################


import numpy as np
import pylab as py

#old code used for different project

#def K10(angle):
#    
#    wave = []
#    
#    wave = np.append(wave,Fresnel(0,1,1.5))
#    wave = np.append(wave,Fresnel(angle,1.5,1))
#    
#    wave = np.append(wave,Fresnel(angle,1,1.5))
#    wave = np.append(wave,Fresnel(angle+Snell(angle,1,1.5),1.5,1))
#    
#    print("flat surface first")
#    print("total s wave")
#    print(wave[0]*wave[2])
#    print("total p wave")
#    print(wave[1]*wave[3])
#    
#    print("angled surface first")
#    print("total s wave")
#    print(wave[4]*wave[6])
#    print("total p wave")
#    print(wave[5]*wave[7])
    
#calculates the S and P polarization transmitance
def Fresnel(theta_i,n1,n2):
    '''
    inputs:
        theta_i = angle of incidence
        n1 = index of refraction of incident material
        n2 = index of refraction of trasmitance material
    returns:
        tuple = [transmission of S polarization, transmission of P polarization]
    '''
    
    theta_i = theta_i *np.pi/180
    
    Rs = (noms(theta_i,n1,n2)/doms(theta_i,n1,n2))* \
         (noms(theta_i,n1,n2)/doms(theta_i,n1,n2))
         
    Rp = (nomp(theta_i,n1,n2)/domp(theta_i,n1,n2))* \
         (nomp(theta_i,n1,n2)/domp(theta_i,n1,n2))
    
    return ([1-Rs,1-Rp])
  
#numerators and denomenators for the Fresnel equations  
def noms(theta_i, n1, n2):
    
    return n1*np.cos(theta_i) - n2*np.sqrt((1-(n1/n2)*np.sin(theta_i))* \
                                            (1-(n1/n2)*np.sin(theta_i)))
    
def doms(theta_i, n1, n2):
    
    return n1*np.cos(theta_i) + n2*np.sqrt((1-(n1/n2)*np.sin(theta_i))* \
                                            (1-(n1/n2)*np.sin(theta_i)))
    
def nomp(theta_i, n1, n2):
    
    return n1*np.sqrt((1-(n1/n2)*np.sin(theta_i))* \
                      (1-(n1/n2)*np.sin(theta_i)))   - n2*np.cos(theta_i)
    
def domp(theta_i, n1, n2):
    
    return n1*np.sqrt((1-(n1/n2)*np.sin(theta_i))* \
                      (1-(n1/n2)*np.sin(theta_i)))   + n2*np.cos(theta_i)

    
#Snells law            
def Snell(theta, n1, n2):
    '''
    inputs:
        theta = angle of incidence
        n1 = index of refraction of initial material
        n2 = index of refraction of final material
    outputs:
        angle of transmitance
    note:
        input angle is in degrees and output angle is in radians
    '''
    
    theta = theta *np.pi/180
    
    return np.arcsin(n1*np.sin(theta)/n2)

#print the lens and calculates ray position on front surface of the lens
def printlens(x_offset, Rcurvature, thickness, diameter):
    '''
    inputs: 
        x_offset = value used to determine y position of ray on the surface of
            front lens
        Rcurvature = curvature of the lens surfaces
        thickness = minimum separation of the lens surfaces
        diameter = length of the lens perpendictular to the optical axis
    outputs:
        plots the lens surfaces
        returns array [y_position of contact between light beam and lens,
                       y_position for calculating angle of incidence,
                       y_position for calculating andle of incidence,
                       array containing y values for back surface]
    '''
    
    front_surface = []
    back_surface = []
    x_cord = [x for x in range(-diameter/2, diameter/2 + 1)]
    
    y_offset = np.sqrt(Rcurvature**2 + (diameter/2)**2)
    
    for i in x_cord:
        
        front_surface.append(gety(i,y_offset,Rcurvature))
        back_surface.append(-gety(i,y_offset,Rcurvature))
    
    front_surface = [x-front_surface[0] for x in front_surface]
    back_surface = [x-back_surface[0] - thickness for x in back_surface]
    
    py.plot(x_cord,front_surface)
    py.plot(x_cord,back_surface)
    py.axis([-diameter/2, diameter/2,-diameter/2, diameter/2])
    
    return([front_surface[diameter/2 + x_offset],
             front_surface[diameter/2 + x_offset - 1],
             front_surface[diameter/2 + x_offset + 1],
             back_surface])

#calculates the y_coordinate of a circle    
def gety(x,y,r):
    
    return -2*y + np.sqrt(4*y**2 - 4*(y**2 + x**2 - r**2))/2
    
#prints the ray and calculates S and P transmittance
def printray(x_offset,Rcurvature, thickness, diameter):
    '''
    inputs:
        x_offset = position of the incoming beam of light
        Rcurvature = curvature of each side of the lens
        thickness = minimum separation of front lens and back lens
        diameter = lenght of the lens
    outputs:
        S and P polarization transmission through lens
        plot of the light beams path
    '''
    
    #call "printlens" to plot the lens surfaces and obtain return values
    y_impact = printlens(x_offset,Rcurvature, thickness, diameter)
    
    #initialize variables
    #y's for input beam
    y_path = [x for x in range(int(np.ceil(y_impact[0])), 100)] 
    #x's for input beam
    x_path = [x_offset for x in y_path] 
    #x's for back surface of lens
    x_check = [x for x in range(-len(y_impact[3])/2,len(y_impact[3])/2+1)] 
    sep = 1000 #storage for separation between data points
    stop = 0 #storage for end point on path of beam
    halt = 0 #storage for end point position on the lens
    n1 = 1 #index of refraction for air
    n2 = 1.5 #index of refraction for glass
    
    #obtain angles for Snells law and use them to obtain slope of light beam
    theta_i = np.arctan((y_impact[1] - y_impact[2])/2)
    theta_t = Snell(theta_i * 180/np.pi, n1, n2)
    slope = 1/np.tan(theta_t)
    
    #for loop plots path of light beam in lens and determines first estimation
    #for contact position of light beam and back surface of the lens
    for x in range(0, 100):
        
        x_path = py.append(x_path, x_offset - x)
        y_path = py.append(y_path, y_impact[0] - slope*x)
        
        for i in range(len(x_check)):
            
            if(np.hypot(x_path[x]-x_check[i-1], y_path[x]-y_impact[3][i-1]) 
               <= sep):
               
               sep = np.hypot(x_path[x]-x_check[i-1],y_path[x]-y_impact[3][i-1])
               stop = x
               halt = i
    
    
    #accurately determins contact point between light beam and back surface lens
    def line(p1, p2):
        A = (p1[1] - p2[1])
        B = (p2[0] - p1[0])
        C = (p1[0]*p2[1] - p2[0]*p1[1])
        return A, B, -C
    
    intersection = True
    i = 0
    while(intersection):
        i += 1
        L1 = line([x_path[stop-i],y_path[stop-i]],
                   [x_path[stop+i],y_path[stop+i]])
        L2 = line([x_check[halt-i],y_impact[3][halt-i]],
              [x_check[halt+i], y_impact[3][halt+i]])
        D  = L1[0] * L2[1] - L1[1] * L2[0]
        Dx = L1[2] * L2[1] - L1[1] * L2[2]
        Dy = L1[0] * L2[2] - L1[2] * L2[0]
        if(D != 0):
            x = Dx / D
            y = Dy / D
            intersection = False
    
    #obtain angles for Snells law and use them to obtain slope of light beam
    #for the back surface of the lens
    theta_f = np.arctan((y_impact[3][halt-i] - y_impact[3][halt+i])/
                         (x_check[halt-i] - x_check[halt+i]))                     
    theta_o = Snell(theta_f * 180/np.pi, n2, n1)
    slope = 1/np.tan(theta_o)
    
    #remove part excess parts of the path arrays
    x_path = x_path[:stop + 1]
    y_path = y_path[:stop + 1]
    
    #for loop creats light beam outside of the back surface of the lens
    for i in range(0, 100):
        
        x_path = py.append(x_path, x - i)
        y_path = py.append(y_path, y - slope*i)
    
    #plot the ray path  
    py.plot(x_path, y_path)                      
    
    #uses Module to determine S and P polarization loss at each surface        
    transmittance = Fresnel(theta_i, n1, n2)
    print "incident s wave transmittance"
    print transmittance[0]
    print "incident p wave transmittance"
    print transmittance[1]
    
    transmittance = Fresnel(theta_f, n2, n1)
    print "transmitted s wave transmittance"
    print transmittance[0]
    print "transmitted p wave transmittance"
    print transmittance[1]            

        
    
    
    
    
    
    
    