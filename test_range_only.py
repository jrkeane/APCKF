from math import pi
import numpy as np
from numpy.ma import zeros
from numpy import genfromtxt, sqrt
import matplotlib.pyplot as plt
import range_only_functions as rof

__author__ = "jim keane"
__date__ = "16 01 20"
__ver__ = "0.1"  # the slow translation prcess
'''
Sanjeev's comments: 
    % File Description: This is the main file for testing range-only tracking
    % Date: 5 Jan 2020
    % Author: Dr Sanjeev Arulampalam
    % DST Group Edinburgh
    % This version simply specifies a quadrant in which the target initial
    % location is in addition to the range-only measurements
'''
'''
Jim's comments
    The goal of this code is to test an algorithm for range-only localisation and homing to a moving target.
    Assumptions:
    Target moving in straight line
    Starting quadrant is known
    Target speed ~1m/s
    Own speed ~2m/s
    Own navigation accuracy is sufficient (e.g. AUV has bottom-lock)
    
    Converting from MATLAB to python by JR Keane, starting 16/01/2020
    
    200 lines converted this afternoon 16 10 2020
    220 to go... 
    12 02 2020 nearly there!


    this is the script from Sanjeevs main test function
    from here we've called the math, as well as enough simulated data to test it (hopefully)

'''

print("APCKF matlab - python conversion")
print("you've got this")
print("Stage 1 - lay out starting parameters "
      "         - set up a config file to manage time, num filters, sampling rate, quadrant, etc")
T = 5 # Sampling time in seconds (time between measurements)
# num_trk = 50
num_trk = 30  # Number of filters used in Angle-Parameterised Cubature Kalman filter (AP-CKF) # how to optimise
t_max = 620  # Duration of scenario in seconds
# N = round(t_max/T)  # Total number of measurements
Xo, Xt = rof.load_ro_scen_data() # load simulation data
N = Xo.shape[1]
# print(N)
# Xo[:,:-1] = []  # why is this here.. not sure its necessary with python, definitely not pythonistic
# print(Xo.shape)
# t = range(N)
# print(t)
sigma_r = 0.10  # standard deviation of range measurements (in meters)
vmax = 3  # Maximum speed of assumed target in m/s
vel_std = vmax/sqrt(3)  # This is the assumed standard deviation of target speed in m/s

q_tild = 0  # This is the process noise intensity parameter in m^2/s^3.

#MC = 100  #Number of Monte Carlo runs
MC = 1  #Number of Monte Carlo runs

rms = np.zeros(N)  # Root Mean Square position error is stored here

print("this example starts in NW quadrant 4 (top-left)")
theta_min = -90*pi/180 # theta_min and theta_max specify the sector in which we expect the initial position of the target
theta_max = 0*pi/180 # these are in radians

print("Now it begins: "
      "iterate through number of monte carlo runs... in this case MC = 1")
for j in range(MC):
    # print(j)
    # if np.mod(j, 2) == 0:  # why are we doing this?
        # print('why we doing this one? ' + str(j))

    print("induce noise into measurements with gen_meas_range_only fucntion")
    Z = rof.gen_meas_range_only(N, sigma_r, Xt, Xo)
    print("")
    print("simulation is set up...")
    print("Call the AP-CKF function to operate on the range-only measurements...")
    print("The output Xthat is a 4 by N matrix containing the target state estimates")
    print("")
    Xt_hat = rof.apckf(Xo,Z,theta_min,theta_max,sigma_r,T,N,vel_std,num_trk,q_tild)

    rof.plot_results(Xt, Xo, Xt_hat)

    # # calc rms position error
    # Xterr = Xthat - Xt
    # x_err = Xterr(0,:)
    # y_err = Xterr(1,:)
    # rms = rms + ((1 / MC) * (np.square(x_err) + np.square(y_err)))

'''
that's the end of the main loop.
sanjeev did some do some extra calcs to compare accuracy against ideal .. 
not sure if this is meant to be separate to the main loop or if it's just been left out as he's only doing one 
monte-carlo iteration so who cares 
'''

# rms = sqrt(rms)
# sig_b = (theta_max - theta_min) / sqrt(12)  # This is the equivalent bearing noise standard deviation
#                                             # for the assumed sector with limits theta_min and theta_max
# sig_r = sigma_r
# sig_vel = vel_std
#
# # The crlb gives the best possible performance for rms position error one
# # could expect from any algorithm given the parameters here
#
# crlb = crlb_range_only(Xt,Xo,sig_b,sig_r,sig_vel,T)

