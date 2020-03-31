from math import pi
import numpy as np
from numpy import genfromtxt, sqrt
import matplotlib.pyplot as plt
from math import pow, ceil, exp

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


# so far translated to line 245 and initialised PfiltM

def load_ro_scen_data():
    '''
    Sanjeev's comments:
        # This loads the target and ownship state vectors
        # target and ownship state vectors are each stored in a 4 by N matrix. Each
        # column represents the state of target (or ownship) with elements (x, y,
        # xdot, ydot) and units (m, m, m/s, m/s)
    '''

    # import where xo = x_own_ship and xt = x_target_ship
    df_xo = genfromtxt('C:/Users/AEA80272/PycharmProjects/gavia/homing/sanjeev_xo.csv', delimiter=',')
    df_xt = genfromtxt('C:/Users/AEA80272/PycharmProjects/gavia/homing/sanjeev_xt.csv', delimiter=',')

    print("import df_xo - ownship data for the simulation:")
    print(df_xo.shape)
    # print("and here is the ownship data: ")
    # print(df_xo)
    print("")
    print("import df_xt - target data for the simulation:")
    print(df_xt.shape)
    # print("and here is the target ship data ")
    # print(df_xt)
    Xo = np.copy(df_xo)
    Xt = np.copy(df_xt)

    return Xo, Xt


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

    r_x = Xt[0, :] - Xo[0, :]
    r_y = Xt[1, :] - Xo[1, :]
    r_diff_sq = np.square(r_x) + np.square(r_y)
    range_true = sqrt(r_diff_sq)  # c^2 = a^2 + b^2, c = sqrt(a^2 + b^2)

    print("array range_true")
    print(range_true)

    range_noise = np.random.rand(1, N)
    print("generate noise")
    range_noise = range_noise * sigma_r
    Z = range_true + range_noise  # add noise to true ranges
    print("generate Z")

    return Z

def apckf(Xo, Z, theta_min, theta_max, sigma_r, T, N, vel_std, num_trk, q_tild):

    '''
    Sanjeevs comments:
        # % file name : apckf.m
        # % file description: This is an Angle-Parameterised Cubature Kalman filter
        # % (AP-CKF) for range-only tracking
        # % Date: 5 Jan 2020
        # % Author: Dr Sanjeev Arulampalam
        # % Organisation: DST Group, Edinburgh
    '''

    print("APCKF function")
    print("Initiliase APCKF "
          "with xo = ownship data, Z = ranges, Q = covariance matrix"
          "wt = weights of APCKF filters")

    R0 = Z[0][0]  # Z gives ranges after noise has been added
    print(R0)  #  RO is first element of R0
    # set up the width of angular region, given num_trk filters
    theta_width = (theta_max - theta_min) / num_trk
    # Now set up theta_vect for the mid-points (in the angular sectors) of all the filters
    theta_vect_min = (theta_min + theta_width / 2)
    theta_vect_max = ((theta_max - theta_width / 2))
    theta_vect = np.arange(theta_vect_min, theta_vect_max, theta_width)
    # init sigb0 as matrix of 1's the size of number of filters to be used
    sigb0 = np.ones(num_trk)
    sigb0 = sigb0 * float(theta_width / sqrt(12))  # the standard deviation in angle for each of the filters
    #  Set up Q the process noise covariance matrix
    Q = np.array([[(1 / 3) * pow(T, 3), 0, (1 / 2) * pow(T, 2), 0],
                  [0, (1 / 3) * pow(T, 3), 0, (1 / 2) * pow(T, 2)],
                  [(1 / 2) * pow(T, 2), 0, T, 0],
                  [0, (1 / 2) * pow(T, 2), 0, T]])  # this matches Sanjeevs

    Q = q_tild * Q  # This is the process noise Covariance matrix
    wt = np.ones(num_trk)    # init wt array
    wt = wt * (1 / num_trk)  # fill with initial weights for the angle-parameterised filters
    x0 = np.zeros((4, num_trk))  # initializing the target states for each of the filters in the AP-CKF
    # print("initializing the target states for each of the filters in the AP-CKF")
    sin_theta_vect = R0 * np.sin(theta_vect)
    cos_theta_vect = R0 * np.cos(theta_vect)
    x0[0, :] = sin_theta_vect + Xo[0][0]
    x0[1, :] = cos_theta_vect + Xo[1][0]
    # print("x0 initial target states for each of the filters. x0 shape and values:")
    # print(x0.shape)
    # print(x0)

    R_var = pow(sigma_r, 2)
    vel_var = pow(vel_std, 2)

    # The following 'for loop' initialises the covariance matrices for each of
    # the num_trk filters
    print("Now run the loop to initialise covariance matrices for each of the num_trk filters")
    first_P0 = True # treat first iter differently
    for i in range(num_trk): # where num_trk is Number of filters used in (AP-CKF)
        # print(i)
        sigma_bet = sigb0[i]
        bet = theta_vect[i]
        sigma_r2 = pow(sigma_r, 2)
        sigma_bet2 = pow(sigma_bet, 2)
        sin_bet2 = pow(np.sin(bet), 2)
        cos_bet2 = pow(np.cos(bet), 2)
        R02 = pow(R0, 2)
        sig_y2 = sigma_r2 * cos_bet2 + R02 * sigma_bet2 * sin_bet2
        sig_x2 = sigma_r2 * sin_bet2 + R02 * sigma_bet2 * cos_bet2
        sig_xy = ((sigma_r2) - (R02) * sigma_bet2) * np.sin(bet) * np.cos(bet)

        PP = np.array([[sig_x2, sig_xy],  # whatever PP stands for... these variable names are killing me
                       [sig_xy, sig_y2]])
        PP1 = np.zeros((2, 2))  # 2 x 2 of zeroes
        PP3 = np.eye(2) * vel_var  # where np.eye gives the identity matrix

        # print("this is the p stack")
        # P_stack = np.dstack((PP,PP1,PP1,PP3))  # this is how we can stack a layer for each i
        # stack the arrays into a 4 x 4. can speed this up later maybe
        top = np.hstack((PP, PP1))
        bottom = np.hstack((PP1, PP3))
        fourfour = np.vstack((top, bottom))
        if first_P0:
            P0 = np.copy(fourfour)
            first_P0 = False

        else:
            P0 = np.dstack((P0, fourfour))

    xf_all = np.copy(x0)
    P_all = np.copy(P0)
    xfm = np.zeros((4))  # xfm will contain the overall state estimate     # xfm translated okay.
    Pfiltm = np.zeros((4, 4))  # Pfiltm the corresponding covariance matrix
    print("Initialised co-variance matrix PFILTM:")
    print(Pfiltm)
    print("")
    print("Form the overall estimate of the AP-CKF based on the Gaussian-Sum formulation of the AP-CKF")

    for i in range(num_trk): # where num_trk is Number of filters used in (AP-CKF)
        # print(i+1)
        wti_x_xfa = wt[i] * xf_all[:, i]
        xfm = np.add(xfm, wti_x_xfa)  # need this xfm to iterate
        # Pfiltm_A = wt[i] * P_all[:, :, i]
        # print(Pfiltm_A)
        # Pfiltm_B = wt[i] * np.matmul(xf_col, xf_row)
        xf_row = np.array(xf_all[:, i], ndmin=2)
        xf_col = xf_row.reshape(-1, 1)
        Pfiltm = Pfiltm + wt[i] * P_all[:, :, i] + wt[i] * np.matmul(xf_col, xf_row)

    xf_transpose = xfm.reshape(-1, 1)  # np.tranpose not working for some arrays
    Pfiltm = Pfiltm - xfm * xf_transpose  # correct to here
    N = Z.shape[1]  # length of Z
    xfs = np.empty((4, N))  # create estimated target states array
    Pfs = np.empty((4, 4, N))  # create corresponding covariance matrices array
    xfs[:, 0] = xfm  # All estimated target states are stored here
    Pfs[:, :, 0] = Pfiltm  # These are the corresponding covariance matrices
    r = pow(sigma_r, 2)  # measurement variance
    A = np.array([[1, 0, T, 0],  # this is the state transition matrix
                  [0, 1, 0, T],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]])

    '''
    Think this is where I need to break out into another function 
    to be able to run real time 
    '''
    print("Run APCKF for rest of measurements, til N. Measurement array shape is:")
    for i in range(1, N):  # where N is total number of measurements available
        wt_sum = 0
        Xo_curr = Xo[:, i]
        Z_curr = Z[0, i]
        for t in range(num_trk): # where num_trk is Number of filters used in (AP-CKF)
            Xthat_curr = xf_all[:, t]
            P_curr = P_all[:, :, t]
            # Now execute the cubature Kalman filter with range-only measurement
            # for the current measurement, given the current state and covariance
            # of the t-th filter
            xtf, Pfilt, likeli = ckf_ro_onestep(Xthat_curr,P_curr,Xo_curr,Z_curr,sigma_r,T,Q)
            xf_all[0,t] = xtf[0]  # not pythonic
            xf_all[1,t] = xtf[1]  # not pythonic
            xf_all[2,t] = xtf[2]  # not pythonic
            xf_all[3,t] = xtf[3]  # not pythonic
            P_all[:,:,t] = Pfilt
            wt[t] = wt[t] * likeli
            wt_sum = wt_sum+wt[t]

        xfm = np.zeros((4))  # xfm will contain the overall state estimate     # xfm translated okay.
        Pfiltm = np.zeros((4, 4))  # Pfiltm the corresponding covariance matrix

        for jj in range(num_trk):
            wt[jj] = wt[jj]/wt_sum
            if wt[jj] < 0.001:
                wt[jj] = 0

            xfm = xfm+wt[jj]*xf_all[:,jj]
            Pfiltm = Pfiltm + wt[jj] * P_all[:,:,jj] + wt[jj] * xf_all[:,jj]*np.transpose(xf_all[:,jj])

        xf_transpose = xfm.reshape(-1, 1)
        Pfiltm = Pfiltm - xfm * xf_transpose
        xfs[:,i] = xfm
        Pfs[:,:,i] = Pfiltm
        # plot_iter(xfs)

    Xt_hat = xfs



    return Xt_hat




def ckf_ro_onestep(Xthat_curr, P_curr, Xo_curr, Z_curr, sigma_r, T, Q):
    '''
    Sanjeevs comments:
        # % This is a Cubature Kalman filter to do range-only tracking
        # % Date: 5 Jan 2020
        # % Author: Dr Sanjeev Arulampalam
        # % Organisation: DST, Australia
    '''

    '''
    Function for CKF_ro_onestep for each iteration
    '''
    r = pow(sigma_r, 2)
    A = np.array([[1, 0, T, 0],  # this is the state transition matrix
                  [0, 1, 0, T],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]])

    nd = 4  # what's ND? size of array?
    two_nd = 2 * nd
    x_pred_pts = np.zeros((nd, two_nd))
    wts = np.ones((two_nd))
    wts = np.multiply((1 / two_nd), wts)
    xp = np.matmul(A, Xthat_curr)
    # xp = xp.reshape(-1,1)
    Pp_A_P = np.matmul(A, P_curr)
    At = np.transpose(A)
    Pp_APAt = np.matmul(Pp_A_P, At)
    Ppred = Pp_APAt + Q
    sqrt_P = np.linalg.cholesky(
        nd * Ppred)  # chol in matlab returns an upper triangular matrix, but this returns the lower as needed

    for j in range(two_nd):
        col_ind = np.mod(j, nd)
        # if col_ind == 0:
        #     col_ind = nd
        j_nd = ceil(j / nd) + 1
        sign_val = pow(-1, j_nd+1)
        # print("sign_val " + str(sign_val))
        check = sqrt_P[:, col_ind]
        x_pred_pts[:, j] = xp + np.multiply(sign_val, check)

    # now compute predicted measurement
    zp_val_vect = np.empty(two_nd)
    for j in range(two_nd):
        v = np.array([x_pred_pts[0, j] - Xo_curr[0], x_pred_pts[1, j] - Xo_curr[1]])
        zp_val = np.linalg.norm(v)
        zp_val_vect[j] = zp_val

    zp = np.matmul(zp_val_vect, np.transpose(wts))

    # compute Pxz and Pzz
    Pxz = np.zeros((nd, 1))
    Pzz = 0

    for j in range(two_nd):
        zp_val = zp_val_vect[j]
        nu = zp_val - zp
        zp_nu = nu
        # wtsjj = wts[j]
        first = np.multiply(wts[j], x_pred_pts[:, j] - xp)
        first = first.reshape(-1, 1)
        zp_nu_t = np.transpose(zp_nu)
        second = np.multiply(first, zp_nu_t)
        Pxz = Pxz + second
        # print(Pxz)
        # pz_f = wts[j] * zp_nu
        Pzz = Pzz + wts[j] * zp_nu * np.transpose(zp_nu)
        # print(Pzz)

    Pzz = Pzz + r
    # he = np.reciprocal(Pzz)
    G = np.multiply(Pxz, np.reciprocal(Pzz))

    nu = Z_curr - zp
    nu_abs = np.abs(nu)

    xp = xp.reshape(-1,1)
    xf = xp + np.multiply(G, nu)
    Gt = np.transpose(G)
    Pfilt = Ppred - G * Pzz * np.transpose(G)

    xtf = xf
    Phat = Pfilt
    low = 1 / sqrt(2 * pi * Pzz)
    up = -0.5 * np.transpose(nu) * np.reciprocal(Pzz) * nu
    likeli = low * exp(up)

    return xtf, Phat, likeli


def plot_results(Xo, Xt, Xt_hat):
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
    target_x = Xt[0, :]
    target_y = Xt[1, :]

    #  Plot ownship track
    print(Xo.shape)
    ownship_x = Xo[0, :]
    ownship_y = Xo[1, :]

    # # plot acpkf estimates
    print(Xt_hat.shape)
    track_x = Xt_hat[0,:]
    track_y = Xt_hat[1,:]

    axs.plot(target_x, target_y, linestyle='solid', label='target', c='r')
    axs.plot(ownship_x, ownship_y, linestyle='solid', label='ownship', c='g')
    axs.plot(track_x,track_y, linestyle='dashed', label='Estimated Track (APCKF)', c='b')
    axs.legend(loc='upper right')
    axs.axis([-600, 400, -600, 400])
    axs.set_aspect('equal')

    plt.show()

    return


def plot_iter(Xt_hat):
    # for figure
    fig, axs = plt.subplots(1, constrained_layout=True)
    fig.canvas.set_window_title('Holthouse')
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized()
    mng.window.activateWindow()
    mng.window.raise_()
    title_str = "Homing - acpkf"
    fig.suptitle(title_str, fontsize="12")

    # #  Plot Target track
    # print(Xt.shape)
    # target_x = Xt[0, :]
    # target_y = Xt[1, :]
    #
    # #  Plot ownship track
    # print(Xo.shape)
    # ownship_x = Xo[0, :]
    # ownship_y = Xo[1, :]

    # # plot acpkf estimates
    print(Xt_hat.shape)
    track_x = Xt_hat[0,:]
    track_y = Xt_hat[1,:]

    # axs.plot(target_x, target_y, linestyle='solid', label='target', c='r')
    # axs.plot(ownship_x, ownship_y, linestyle='solid', label='ownship', c='g')
    axs.plot(track_x,track_y, linestyle='dashed', label='Estimated Track (APCKF)', c='b')
    axs.legend(loc='upper right')
    axs.axis([-400, 200, -400, 200])
    axs.set_aspect('equal')

    plt.show()

    return

def crlb_range_only(Xt, Xo, sig_b, sig_r, sig_vel, T):
    return crlb
