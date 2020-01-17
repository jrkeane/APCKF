from math import pi
import numpy as np
from numpy import genfromtxt, sqrt
import matplotlib.pyplot as plt
from math import pow

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
    # print(df_xo)
    print("")
    print("import df_xt - target data for the simulation")
    print(df_xt.shape)
    print("and here is the target ship data ")
    # print(df_xt)

    Xo = np.copy(df_xo)
    # print(Xo)
    Xt = np.copy(df_xt)

    return Xo, Xt

def apckf(Xo,Z,theta_min,theta_max,sigma_r,T,N,vel_std,num_trk,q_tild):
    #     % file name : apckf.m
    # % file description: This is an Angle-Parameterised Cubature Kalman filter
    # % (AP-CKF) for range-only tracking
    # % Date: 5 Jan 2020
    # % Author: Dr Sanjeev Arulampalam
    # % Organisation: DST Group, Edinburgh

    '''initialise '''

    print("heres Z (array of ranges)")
    print(Z)
    print("get first element from Z - just to check ")
    R0 = Z[0][0]  # this is slightly different to sanjeevs cause randn
    print(R0)
    theta_width = (theta_max-theta_min)/num_trk  # the width of angular region, given num_trk filters
    print("theta width (width of angular region, given filters ")
    print(theta_width)

    theta_vect_min = (theta_min+theta_width/2)
    theta_vect_max = ((theta_max-theta_width/2))
    theta_vect = np.arange(theta_vect_min, theta_vect_max, theta_width)  # the mid-points (in the angular sectors)
                                                                                # of all the filters
    # print("theta vector")
    # print(theta_vect)
    sigb0 = np.ones(num_trk)
    sigb0 = sigb0 * float(theta_width/sqrt(12))  # the standard deviation in angle for each of the filters
    # print("sigb0")
    # print(sigb0)

    Q = np.array([[(1/3)*pow(T,3), 0, (1/2)*pow(T,2), 0],
                  [0, (1/3)*pow(T,3), 0, (1/2)*pow(T,2)],
                  [(1/2)*pow(T,2), 0, T, 0],
                  [0, (1/2)*pow(T,2), 0, T]])  # this matches Sanjeevs

    Q = q_tild*Q  #This is the process noise Covariance matrix
    # print(Q)# this matches Sanjeevs


    wt = np.ones(num_trk)
    wt = wt * (1/num_trk)  # These are the initial weights for the angle-parameterised filters
    # print("These are the initial weights for the angle-parameterised filters")
    # print(wt)
    x0 = np.zeros((4,num_trk))  # initializing the target states for each of the filters in the AP-CKF
    print("initializing the target states for each of the filters in the AP-CKF")
    # print(x0)
    # print()
    print("")
    sin_theta_vect = R0*np.sin(theta_vect)
    # print(sin_theta_vect.shape)
    # print(sin_theta_vect)
    cos_theta_vect = R0*np.cos(theta_vect)
    # print("Xown")
    # print(Xo[0][0])
    # print(x0[1,:])
    x0[0,:] = sin_theta_vect+Xo[0][0]
    x0[1,:] = cos_theta_vect+Xo[1][0]
    # x0[3:4,:] = [0,0]*np.ones(1,num_trk)
    print("x0")
    print(x0.shape)
    print(x0)


    R_var = pow(sigma_r,2)
    vel_var = pow(vel_std,2)
    '''
    translated up to here 17 01 20 
    '''
    # The following 'for loop' initialises the covariance matrices for each of
    # the num_trk filters
    print("Now run the loop to initialise covariance matrices for each of the num_trk filters")
    first_P0 = True


    for i in range(num_trk):
        print(i)
        sigma_bet =sigb0[i]
        bet = theta_vect[i]
        sigma_r2 = pow(sigma_r,2)
        sigma_bet2 = pow(sigma_bet,2)
        sin_bet2 = pow(np.sin(bet),2)
        cos_bet2 = pow(np.cos(bet),2)
        R02 = pow(R0,2)
        sig_y2 = sigma_r2*cos_bet2+R02*sigma_bet2*sin_bet2
        sig_x2 = sigma_r2*sin_bet2+R02*sigma_bet2*cos_bet2
        sig_xy = ((sigma_r2)-(R02)*sigma_bet2)*np.sin(bet)*np.cos(bet)

        # we want a 4 x 4 x num_trk matrix  i.e np.array([PP, np.zeros(2,2)],
        #                                               [np.zeros(2,2), vel_var*np.eye(2)]])

        PP = np.array([[sig_x2, sig_xy],  # whatever PP stands for... these variables are killing me
                        [sig_xy, sig_y2]])
        PP1 = np.zeros((2,2))    # 2 x 2 of zeroes
        PP3 = np.eye(2)*vel_var  # where np.eye gives the identity matrix

        # print("Now let's figure out these matrices")
        # print("this is the p stack")
        # P_stack = np.dstack((PP,PP1,PP1,PP3))  # this is how we can stack a layer for each i
        # print(P_stack.shape)
        # print(P_stack)
        print("")
        # stack the arrays into a 4 x 4. can speed this up heaps later
        top = np.hstack((PP,PP1))
        bottom = np.hstack((PP1,PP3))
        # print("top")
        # print(top)
        fourfour = np.vstack((top,bottom))
        if first_P0:
            P0 = np.copy(fourfour)
            print("this is the four x four stack")
            print(P0.shape)
            print(P0)
            first_P0 = False
        else:
            P0 = np.dstack((P0, fourfour))
            print(P0.shape)

    print(P0)

    xf_all = np.copy(x0)
    P_all = np.copy(P0)


'''
% Now form the overall estimate of the AP-CKF based on the Gaussian-Sum
% formulation of the AP-CKF
% xfm will contain the overall state estimate and Pfiltm the corresponding
% covariance matrix
'''

xfm = np.zeros((4,1))
Pfiltm = np.zeros(4,4)
for i in range(num_trk):
    xfm = xfm + wt[i] * xf_all[:,i]
    Pfiltm = Pfiltm+wt[i] * P_all(:,:,i) + wt[i] * xf_all[:,i] * xf_all(:,i)';

Pfiltm=Pfiltm-xfm*xfm';

N = length(Z);

xfs = nan(4,N);
Pfs = nan(4,4,N);

''' Now run the AP-CKF for the rest of the measurements (2 till N) '''

    Xt_hat = 0

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

    range_noise = np.random.rand(1,N)
    print("generate noise")
    # print(range_noise)
    range_noise = range_noise * sigma_r
    # print(range_noise)
    # print("hold")
    Z = range_true + range_noise  # add noise to true ranges
    print("generate Z")
    # print(Z)

    return Z

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

