from math import pi
import numpy as np
from numpy.ma import zeros
from numpy import genfromtxt, sqrt
import matplotlib.pyplot as plt

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

'''


def load_ro_scen_data():
    '''
    Sanjeev's comments:
    # This loads the target and ownship state vectors
    # target and ownship state vectors are each stored in a 4 by N matrix. Each
    # column represents the state of target (or ownship) with elements (x, y,
    # xdot, ydot) and units (m, m, m/s, m/s)
    '''

    # let's try read straight in with pandas
    # df_xo = pd.read_excel(r'C:\Users\AEA80272\PycharmProjects\gavia\homing\sanjeev_xo.xlsx')
    # df_xt = pd.read_excel(r'C:\Users\AEA80272\PycharmProjects\gavia\homing\sanjeev_xt.xlsx')

    # cancel that,
    # import where xo = x_own_ship and xt = x_target_ship
    df_xo = genfromtxt('C:/Users/AEA80272/PycharmProjects/gavia/homing/sanjeev_xo.csv', delimiter=',')
    df_xt = genfromtxt('C:/Users/AEA80272/PycharmProjects/gavia/homing/sanjeev_xt.csv', delimiter=',')

    print("import df_xo - ownship data for the simulation")
    print(df_xo.shape)
    print("and here is the ownship data ")
    print(df_xo)
    print("")
    print("import df_xt - target data for the simulation")
    print(df_xt.shape)
    print("and here is the target ship data ")
    print(df_xt)

    Xo = np.copy(df_xo)
    Xt = np.copy(df_xt)

    return Xo, Xt

def apckf():

    return Xt_hat

def gen_meas_range_only(N, sigma_r, Xt, Xo):
    '''
    Sanjeev's comments:
    file description: This generates the range - only measurements
    Given Xo, and Xt, generates a vector Z of measurements;
    Author: Dr Sanjeev Arulampalam
    Organisation: DST Group Edinburgh
    Date: 5 Jan 2020
    '''

    Z = np.zeros(N)  # initialise array Z for range only (with noise)

    # xrel = Xt - Xo  # get the difference between target state and ownship state
    r_x = Xt[0,:] - Xo[0,:]
    # print(r_x)
    # print(r_x.shape)
    r_y = Xt[1,:] - Xo[1,:]
    # print(r_y)

    r_diff_sq = np.square(r_x) + np.square(r_y)
    range_true = sqrt(r_diff_sq)  # c^2 = a^2 + b^2, c = sqrt(a^2 + b^2)

    print("array range_true")
    print(range_true)

    Z = range_true + (sigma_r * np.random.randn(0, N))  # add noise to true ranges

    return Z


def apckf(Xo,Z,theta_min,theta_max,sigma_r,T,N,vel_std,num_trk,q_tild):

    return Xthat

def ckf_ro_onestep():

    return

def plot_results(Xo, Xt):
    # for figure
    fig, axs = plt.subplots(1, constrained_layout=True)
    fig.canvas.set_window_title('Holthouse')
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized()
    mng.window.activateWindow()
    mng.window.raise_()
    title_str = "Homing - acpkf"
    fig.suptitle(title_str, fontsize="12")

    #  Plot Target track
    print(Xt.shape)
    target_x = Xt[0,:]
    target_y = Xt[1,:]

    #  Plot ownship track
    print(Xo.shape)
    ownship_x = Xo[0,:]
    ownship_y = Xo[1,:]

    # # plot acpkf estimates
    # print(Xt_hat.shape)
    # track_x = Xt_hat[0,:]
    # track_y = Xt_hat[1,:]

    axs.plot(target_x, target_y,  linestyle='solid', label='target', c='r')
    axs.plot(ownship_x, ownship_y, linestyle='solid', label='ownship', c='g')
    # axs.plot(track_x,track_y, linestyle='solid', label='Estimated Track (APCKF)', c='b')
    axs.legend(loc='upper right')
    axs.axis([-400, 200, -400, 200])
    axs.set_aspect('equal')

    plt.show()

    return


def crlb_range_only(Xt,Xo,sig_b,sig_r,sig_vel,T):

    return crlb


'''

this is the script from Sanjeevs main test function
from here we've called the math, as well as enough simulated data to test it (hopefully)

'''

T = 5 # Sampling time in seconds (time between measurements)
# num_trk = 50
num_trk = 30  # Number of filters used in Angle-Parameterised Cubature Kalman filter (AP-CKF) # how to optimise
t_max = 620  # Duration of scenario in seconds
N = round(t_max/T)  # Total number of measurements
print(N)

Xo, Xt = load_ro_scen_data() # load simulation data
N = Xo.shape[1]
print(N)

# Xo[:,:-1] = []  # why is this here.. not sure its necessary with python, definitely not pythonistic
print(Xo.shape)

# t = range(N)
# print(t)

sigma_r = 10  # standard deviation of range measurements (in meters)

# vel_std = 1.5  # This is the assumed standard deviation of target speed in m/s
vmax = 3  # Maximum speed of assumed target in m/s
vel_std = vmax/sqrt(3)  # This is the assumed standard deviation of target speed in m/s

q_tild = 0  # This is the process noise intensity parameter in m^2/s^3.

#MC = 100  #Number of Monte Carlo runs
MC = 1  #Number of Monte Carlo runs

rms = np.zeros(N)  # Root Mean Square position error is stored here

# specify quadrant - in this case Quad 4 (top-left):
theta_min = -90*pi/180 # theta_min and theta_max specify the sector in which we expect the initial position of the target
theta_max = 0*pi/180 # these are in radians

for j in range(MC):
    print(j)
    if np.mod(j, 2) == 0:  # why are we doing this?
        print('why we doing this one? ' + str(j))

    Z = gen_meas_range_only(N, sigma_r, Xt, Xo)

    '''
    # Now call the AP-CKF function to operate on the range-only measurements
    # The output Xthat is a 4 by N matrix containing the target state estimates
    '''

    # Xt_hat = apckf(Xo,Z,theta_min,theta_max,sigma_r,T,N,vel_std,num_trk,q_tild)

    plot_results(Xt, Xo)

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

