import numpy as np
import matplotlib.pyplot as plt
import re
import argparse
from scipy import fft 
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
import subprocess
import statistics as st
import sys
import os
import csv
# Parse Arguments from terminal (README)
# ex: >> python3 LR_script.py data_Wed_Jul_10_15_22_53_2024 twiss.out_pp24-100GeV-e1store Blue --no-NLPlot
def argparser(args):
    parser = argparse.ArgumentParser(description='Use Phase Transfer Matrix to find Optics at Drift Spaces')
    parser.add_argument('sdds_file', type=str,
                        help='path to RHIC sdds file (ex: /operations/app_store/RunData/[run year]/[run number]/RHIC/Orbit/TBT/[ring]/[sdds file]) or name of converted sdds file (name of txt file in your directory)')
    parser.add_argument('twiss_file', type=str,
                        help='path to RHIC twiss file for optics comparison; ex: /operations/app_store/Ramps/MADx/[beam]/[twiss file]')
    parser.add_argument('ring', type=str,
                        help='Color of Ring used (Blue or Yellow)')
    # parser.add_argument('--isPSPlot', type=bool, required = False, help='Plot Phase Space graphs; Default = False')
    parser.add_argument('--no-NLPlot', dest='isNLPlot',
                        action='store_false',
                        help='Disable Nonlinearity graphs; Default is to enable')

    parser.set_defaults(isNLPlot=True)

    return parser.parse_args(args)

# Extract data from converted sdds file
# Data Output: BPM Values, Names of bpms, BPM Positions 
def data_from_sdds(sdds_file, ring):
    Name_xbpms = []
    Name_ybpms = []

    # if it is an sdds file, convert to ascii
    if sdds_file[-5:] == ".sdds":
        index_ring = sdds_file.find(ring)
        index_sdds = sdds_file.find('.sdds')
        sdds_file = "data_" + sdds_file[index_ring: index_sdds].replace(":", "_")
        if not os.path.exists(sdds_file):
            subprocess.run(["sddsconvert", "-ascii", sdds_file, sdds_file])

    with open("TBT_Data/" + sdds_file) as output:
        Lines = output.readlines()
        index_ix = index_jx = 0
        index_iy = index_jy = 0
        
        # retrieve # BPMs and turn number
        N_BPMx = int(Lines[63])
        N_BPMy = int(Lines[65])
        N_turns = int(Lines[64])
        Pos_xbpms = np.zeros([N_turns, N_BPMx])
        Pos_ybpms = np.zeros([N_turns, N_BPMy])
        S_xbpms = np.zeros([N_turns, N_BPMx])
        S_ybpms = np.zeros([N_turns, N_BPMy])
        
        init_i = 414
        end_y_index = N_BPMy*N_turns + init_i + 1
        for i, line in enumerate(Lines[1:]):
            Line = []
            #xBPM Names
            if i >= 67 and i < 95:
                for word in re.split("\s+", line):
                    Line.append(word)
                Name_xbpms.append(Line[:-1])
            
            #yBPM Names
            if i >= 241 and i < 269:
                for word in re.split("\s+", line):
                    Line.append(word)
                Name_ybpms.append(Line[:-1])
            
            #s and positions for xBPM
            if i > init_i:
                row = re.split("\s+", line)
                S_xbpms[index_ix][index_jx] = float(row[1])
                Pos_xbpms[index_ix][index_jx] = float(row[2])
                index_jx = (i - init_i)%N_BPMx

                # if i < init_i + N_BPMx and float(row[5]) != 0:
                #     Bad_xBpms.append(index_jx - 1)
                
                if (i - init_i)%N_BPMx == 0 and (i - init_i) != 0:
                    index_ix += 1
            
            #s and positions for yBPM
            if i > init_i and i < end_y_index:
                row = re.split("\s+", line)
                S_ybpms[index_iy][index_jy] = float(row[6])
                Pos_ybpms[index_iy][index_jy] = float(row[7])
                index_jy = (i - init_i)%N_BPMy

                # if i < init_i + N_BPMy and float(row[10]) != 0:
                #     Bad_yBpms.append(index_jy - 1)
                
                if (i - init_i)%N_BPMy == 0 and (i - init_i) != 0:
                    index_iy += 1
            
    Name_xbpms = [item for row in Name_xbpms for item in row]
    Name_ybpms = [item for row in Name_ybpms for item in row]
    # use only 1024 turns
    N_turns = 1024
    BPMx = Pos_xbpms[:N_turns].T
    BPMy = Pos_ybpms[:N_turns].T
    return BPMx, BPMy, Name_xbpms, Name_ybpms, S_xbpms, S_ybpms#, Bad_xBpms, Bad_yBpms

# Extract data from twiss file
# Data Output: Model Beta and Model Alpha Values
def data_from_twiss(twiss_file):
    twiss_file = "Twiss_files/" + twiss_file
    BetxModel = {}
    AlfxModel = {}
    BetyModel = {}
    AlfyModel = {}
    sModel = {}
    # Get Linear optics from MADx
    with open(twiss_file) as output:
        Lines = output.readlines()
        for line in Lines:
            linelist = re.split("\s+", line)
            # print(linelist)
            if linelist[0] == '*':
                idx_s = linelist.index('S')
                idx_BETX = linelist.index('BETX')
                idx_ALFX = linelist.index('ALFX')
                idx_BETY = linelist.index('BETY')
                idx_ALFY = linelist.index('ALFY')
            
            if linelist[0] == '':
                sModel[linelist[1][1:-1]] = float(linelist[idx_s])
                BetxModel[linelist[1][1:-1]] = float(linelist[idx_BETX])
                AlfxModel[linelist[1][1:-1]] = float(linelist[idx_ALFX])
                BetyModel[linelist[1][1:-1]] = float(linelist[idx_BETY])
                AlfyModel[linelist[1][1:-1]] = float(linelist[idx_ALFY])
    BetxBPMModel = []
    BetyBPMModel = []
    AlfxBPMModel = []
    AlfyBPMModel = []
    # Get MADx values of Beta and Alpha from BPMs
    for elem in sModel:
        if ("_B" in elem and 
            elem != "G10_BX.1" and 
            elem != "G10_BX.2" and 
            "Q" not in elem and 
            "T" not in elem and 
            "BB" not in elem):
            if "V" not in elem: 
                BetxBPMModel.append(BetxModel[elem])
                AlfxBPMModel.append(AlfxModel[elem])
            if "H" not in elem and elem != "BI9_B3.1": 
                BetyBPMModel.append(BetyModel[elem])
                AlfyBPMModel.append(AlfyModel[elem])
    return BetxBPMModel, BetyBPMModel, AlfxBPMModel, AlfyBPMModel, BetxModel, AlfxModel, BetyModel, AlfyModel

# Plots TBT data and Tune for first BPM
# Output: Tunes, user inputs for LR analysis
def initial_plots(BPMx, BPMy):
    # plt.ion()  # Turn on interactive mode
    # Look at first BPM Signals
    BPMx0 = BPMx[0]
    N_turns = len(BPMx0)

    BPMx0_mean = np.mean(BPMx0)
    BPMx0_std = np.std(BPMx0)
    print()
    print("Mean of First Horizontal BPM =", BPMx0_mean)
    print("Standard deviation of First Horizontal BPM =", BPMx0_std)

    BPMy0 = BPMy[0]
    BPMy0_mean = np.mean(BPMy0)
    BPMy0_std = np.std(BPMy0)
    print("Mean of First Vertical BPM =", BPMy0_mean)
    print("Standard deviation of First Vertical BPM =", BPMy0_std)

    fig, Ax = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    x = np.linspace(0, N_turns, N_turns)
    Ax[0][0].plot(x, BPMx0)
    Ax[0][0].axhline(y = BPMx0_mean, color = 'b', linestyle = '-', linewidth = .75, label = r"mean = %f"%BPMx0_mean)
    Ax[0][0].set_xlabel("Turns")
    Ax[0][0].set_ylabel("Horizontal Offset [mm]")
    Ax[0][0].set_title("Horizontal Signal of BPM 1")
    Ax[0][0].legend()

    Ax[0][1].plot(x, BPMy0)
    Ax[0][1].axhline(y = BPMy0_mean, color = 'b', linestyle = '-', linewidth = .75, label = r"mean = %f"%BPMy0_mean)
    Ax[0][1].set_xlabel("Turns")
    Ax[0][1].set_ylabel("Vertical Offset [mm]")
    Ax[0][1].set_title("Vertical Signal of BPM 1")
    Ax[0][1].legend()

    # Look at FFT graphs and tune of beam
    Omega = np.linspace(0, N_turns//2, N_turns//2)/N_turns
    BPMX0 = np.abs(fft.rfft(BPMx0)[1:])
    BPMx0_Tune_scipy = Omega[np.argmax(BPMX0)]
    Ax[1][0].plot(Omega, BPMX0)
    Ax[1][0].axvline(x = BPMx0_Tune_scipy, color = 'b', linestyle = '-', linewidth = .75, label = r"$\nu_x$ = %f"%BPMx0_Tune_scipy)
    Ax[1][0].set_title("FFT at Horizontal BPM 1")
    Ax[1][0].set_xlabel(r"$\omega$")
    Ax[1][0].set_ylabel(r"F($\omega$)")

    BPMY0 = np.abs(fft.rfft(BPMy0)[1:])
    BPMy0_Tune_scipy = Omega[np.argmax(BPMY0)]
    Ax[1][1].plot(Omega, BPMY0)
    Ax[1][1].axvline(x = BPMy0_Tune_scipy, color = 'b', linestyle = '-', linewidth = .75, label = r"$\nu_y$ = %f"%BPMy0_Tune_scipy)
    Ax[1][1].set_title("FFT at Vertical BPM 1")
    Ax[1][1].set_xlabel(r"$\omega$")
    Ax[1][1].set_ylabel(r"F($\omega$)")
    # plt.show(block=False)
    plt.show()
    # print()

    # Based on the plots, determine parameters for LR
    #  Provide valid starting turn and number of turns used in LR
    #  Event loop to allow for non-blocking interaction
    # while True:
    #     try:
    #         start_turn_x = input("Starting turn for Horizontal (default 25): ")
    #         start_turn_x = 25 if start_turn_x == '' else int(start_turn_x)
    #         start_turn_y = input("Starting turn for Vertical (default 525): ")
    #         start_turn_y = 525 if start_turn_y == '' else int(start_turn_y)

    #         print("Vary # turns while keeping starting turn constant:")
    #         start_interval = input("Number of turns to start for Linear Regression (default 50): ")
    #         start_interval = 50 if start_interval == '' else int(start_interval)
    #         end_interval = input("Number of turns to end for Linear Regression (default 200): ")
    #         end_interval = 200 if end_interval == '' else int(end_interval)

    #         # print("Vary starting turns while keeping # turns constant (default 100): ")
    #         # const_interval = 100 if const_interval == '' else int(const_interval)
    #         plt.ioff()  # Turn off interactive mode
    #         break  # Break the loop once inputs are provided
    #     except ValueError as e:
    #         print("Invalid input, please enter numeric values.")
    #     except Exception as e:
    #         print(f"An unexpected error occurred: {str(e)}")
    #         break

    # plt.close(fig)  # Close the figure
    # start_interval //= 10
    # end_interval //= 10
    start_turn_x, start_turn_y, start_interval, end_interval = 25, 525, 5, 20
    return start_turn_x, start_turn_y, start_interval, end_interval#, const_interval

# Look at all FFTs to figure out which BPMs were bad; ones in sdds file sometimes are not correct
# Output: Bad bpms in horizontal and vertical direction
def Tune_Map(BPMx, BPMy):
    N_BPMx = len(BPMx)
    N_BPMy = len(BPMy)
    N_turns = len(BPMx[0])
    Omega = np.linspace(0, N_turns//2, N_turns//2)/N_turns

    BPMx_Tune_Map = np.zeros(N_BPMx)
    for i in range(N_BPMx):
        BPMX = np.abs(fft.rfft(BPMx[i])[1:]) 
        BPMx_Tune_Map[i] = Omega[np.argmax(BPMX)]
    
    BPMy_Tune_Map = np.zeros(N_BPMy)
    for i in range(N_BPMy):
        BPMY = np.abs(fft.rfft(BPMy[i])[1:])
        BPMy_Tune_Map[i] = Omega[np.argmax(BPMY)]

    nu_x, nu_y = st.mode(BPMx_Tune_Map), st.mode(BPMy_Tune_Map)
    print("Average Horizontal tune:", nu_x)
    print("Average Vertical tune:", nu_x)

    plt.plot(BPMx_Tune_Map, color = 'b', label = r"$\nu_x$")
    plt.axhline(y = nu_x, color = 'b', linestyle = '-', linewidth = .75, label = r"$\nu_x$ = %f"%nu_x)
    plt.plot(BPMy_Tune_Map, color = 'r', label = r"$\nu_y$")
    plt.axhline(y = nu_y, color = 'r', linestyle = '-', linewidth = .75, label = r"$\nu_y$ = %f"%nu_y)
    plt.xlabel("BPM number")
    plt.ylabel(r"$\nu$")
    plt.title("Tune measurements of BPMS")
    plt.legend()
    plt.show()

    def Bad_Bpms(Tune_Map, nu, N_BPM):
        Tune_Offset = np.abs(Tune_Map - np.ones(N_BPM)*nu)
        Bad_BPM = []
        for i, offset in enumerate(Tune_Offset):
            if offset > .005:
                Bad_BPM.append(i)
        return Bad_BPM
    Bad_xBPMs = Bad_Bpms(BPMx_Tune_Map, nu_x, N_BPMx)
    Bad_yBPMs = Bad_Bpms(BPMy_Tune_Map, nu_y, N_BPMy)
    return nu_x, nu_y, Bad_xBPMs, Bad_yBPMs

# Does Linear Regression calculation on TBT data to calculate Phase Transfer Matrix
# Data Output: Twiss parameters and beta beat
def LR_calculation(Interval, Bpms, BPM, S_bpms, Name, BetBPMModel, tune_scipy, ring, unExp = True):
    bpm1, bpm2 = Bpms
    b1_m, b2_m = BetBPMModel[bpm1], BetBPMModel[bpm2]
    Data1 = [Name[bpm1], bpm1]
    Data2 = [Name[bpm2], bpm2]
    
    # Transforms tbt data from exp to pure sinusoidal
    def transformtbt(tbt):
        n_turns = len(tbt)
        n_more = n_turns + 1
        x = np.linspace(0, interval_turn - 1, interval_turn)
        # Define the exponential function to fit
        b = -0.006573250198
        def exponential_func(x, a, b, c, d):
            return a * np.exp(b * x)*np.cos(2*np.pi*c *x + d)
        Amp = (np.max(tbt) - np.min(tbt))/2
        popt, pcov = curve_fit(exponential_func, x, tbt, [Amp, b, tune_scipy, 0])
        a_fit, b_fit, c_fit, d_fit = popt
    
        def exponential_func_mod(x, a, c, d):
            return a * np.cos(2*np.pi*c *x + d)
        x_more = np.linspace(0, n_more - 1, n_more)
        fitted_curve_mod = exponential_func_mod(x_more, a_fit, c_fit, d_fit)
        return fitted_curve_mod
    
    #Organize position and speed data
    initial_turn, interval_turn = Interval
        
    # Scaled position and speed between two bpms
    C = 3833.845181
    if S_bpms[0, bpm1] > S_bpms[0, bpm2]: 
        # Special treatment for IP6:
        L = S_bpms[0, bpm2] - S_bpms[0, bpm1] + C
        if ring == "Blue":
            x_data1_all = BPM[bpm1][:-1] - np.mean(BPM[bpm1])
            x_data2_all = BPM[bpm2][1:] - np.mean(BPM[bpm2])
        elif ring == "Yellow":
            x_data1_all = BPM[bpm1][1:] - np.mean(BPM[bpm1])
            x_data2_all = BPM[bpm2][:-1] - np.mean(BPM[bpm2])
        else:
            print("invalid ring (Blue or Yellow)")
            sys.exit()
    else: 
        L = S_bpms[0, bpm2] - S_bpms[0, bpm1]
        x_data1_all = BPM[bpm1] - np.mean(BPM[bpm1])
        x_data2_all = BPM[bpm2] - np.mean(BPM[bpm2])
    xp_data_all = (x_data2_all - x_data1_all)/L

    # variable data: x_data (bpm, var) (current turn)

    n_turns = 450
    x_data11 = x_data1_all[initial_turn: initial_turn + interval_turn]
    x_data12 = xp_data_all[initial_turn: initial_turn + interval_turn]
    
    x_data21 = x_data2_all[initial_turn: initial_turn + interval_turn]

    if unExp:
        x_data11_trans = transformtbt(x_data11)
        x_data12_trans = transformtbt(x_data12)
        x_data21_trans = transformtbt(x_data21)

    x_data11 = x_data11_trans[:interval_turn]
    x_data12 = x_data12_trans[:interval_turn]
    x_data21 = x_data21_trans[:interval_turn]
    x_data22 = np.copy(x_data12)
    
    X_data1 = np.vstack((x_data11, x_data12)).T
    X_data2 = np.vstack((x_data21, x_data22)).T

    # Y Train data (next turn)
    y_data11 = x_data11_trans[1:interval_turn + 1]
    y_data12 = x_data12_trans[1:interval_turn + 1]

    y_data21 = x_data21_trans[1:interval_turn + 1]
    y_data22 = np.copy(y_data12)

    Y_data1 = np.vstack((y_data11, y_data12)).T
    Y_data2 = np.vstack((y_data21, y_data22)).T
            
    #LR Calculation:
    # Y = MX; M: 2x2 (doesn't account for coupling)
    result1 = LinearRegression().fit(X_data1, Y_data1)
    result2 = LinearRegression().fit(X_data2, Y_data2)
    
    # Obtain Phase Transfer matrix from LR
    def Mat_param(result):
        M = result.coef_
        det = np.linalg.det(M)
        # Correct for the determinant to be symplectic
        M /= (np.sqrt(det))
        return M
    
    M1 = Mat_param(result1)
    M2 = Mat_param(result2)

    def Twiss_calc(M):
        if M[0, 1] > 0:
            phi = np.arccos((M[0, 0] + M[1, 1])/2)
        else: phi = 2*np.pi - np.arccos((M[0, 0] + M[1, 1])/2)
        beta = M[0, 1]/np.sin(phi)
        alpha = (M[0, 0] - M[1, 1])/(2*np.sin(phi))
        tune = phi/(2*np.pi)
        Twiss = np.array([tune, phi, alpha, beta])
        return Twiss
    
    # Obtain Twiss Parameters from Phase Transfer Matrix
    # print(Name[bpm1], Name[bpm2], bpm1, bpm2, b1_m, b2_m)
    Twiss1 = Twiss_calc(M1)
    Twiss2 = Twiss_calc(M2)
    tune1, phi1, alf1, bet1 = Twiss1
    tune2, phi2, alf2, bet2 = Twiss2
    # print(M1)
    Data1.extend(Twiss1)
    Data2.extend(Twiss2)
    
    # beta beat at bpm 1, 2:
    beat1 = (b1_m - bet1)/b1_m
    beat2 = (b2_m - bet2)/b2_m
    Data1.extend([b1_m, beat1])
    Data2.extend([b2_m, beat2])

    return Data1, Data2

# Calculates Optics at all drifts
def Compute_Optics(TurnList, Drift_bpm, Name_bpms, BetBPMModel, BPM, S_bpms, tune, ring):
    #  Collect data using LR
    Data = np.zeros([len(TurnList), len(Drift_bpm)*2, 7])
    BPMDriftNames = []
    for i, interval in enumerate(TurnList):
        print(str(interval[0]), " to ", str(interval[1]), " turns")

        unexp = True
        for j in range(len(Drift_bpm)):
            Bpms = Drift_bpm[j]
            Interval = [interval[0], interval[1] - interval[0]]

            # You need both BPMs to get twiss parameters, if even one BPM is bad, both are bad.
            # print(Bpms)
            if None in Bpms:
                bpm1, bpm2 = Drift_bpm[j]
                Data1 = [Name_bpms[bpm1], bpm1, 0, 0, 0, 0, BetBPMModel[bpm1], 1]
                Data2 = [Name_bpms[bpm2], bpm2, 0, 0, 0, 0, BetBPMModel[bpm2], 1]
            else:
                Data1, Data2 = LR_calculation(Interval, Bpms, BPM, S_bpms, Name_bpms, BetBPMModel, tune, ring, unexp)

            if i == 0: 
                BPMDriftNames.append(Data1[0])
                BPMDriftNames.append(Data2[0])
            Data[i][j*2] = np.array(Data1[1:])
            Data[i][j*2 + 1] = np.array(Data2[1:])

    # Obtain relevant information
    Beat_Matrix = Data[:, :, -1]
    Bet_Matrix = Data[:, :, -3]
    Alf_Matrix = Data[:, :, -4]

    return BPMDriftNames, Beat_Matrix, Bet_Matrix, Alf_Matrix

# Plots Optics with respect to # of turns used in LR
# Output: mean and standard deviation of optics
def NonlinearityGraph(M1, M2, BPMNames, title, isNLPlot, Model):
    assert(len(M1[0]) == len(M2[0]))
    nrow = 6
    ncol = 4
    if isNLPlot:
        fig, Ax = plt.subplots(nrow, ncol*2, gridspec_kw = {'width_ratios':[4, 1]*ncol}, figsize=(18, 10), layout='constrained', dpi = 75)
        fig.suptitle(title)
    
    Center1 = []
    Stddev1 = []
    Center2 = []
    Stddev2 = []
    N = len(M1[0])
    X_axis = np.linspace(1, N, N)
    
    for i in range(nrow):
        for j in range(ncol):
            Matrix_target1 = M1[ncol*i + j]
            center1, stddev1 = np.mean(Matrix_target1), np.std(Matrix_target1)
            Center1.append(center1)
            Stddev1.append(stddev1)

            Matrix_target2 = M2[ncol*i + j]
            center2, stddev2 = np.mean(Matrix_target2), np.std(Matrix_target2)
            Center2.append(center2)
            Stddev2.append(stddev2)

            Avg_val = (center1 + center2)/2
            
            if isNLPlot:
                # Plot optics 
                Ax[i][j*2].set_title(BPMNames[ncol*i + j])
                Ax[i][j*2].plot(X_axis, Matrix_target1, '-bo', ms = 3)
                Ax[i][j*2].axhline(y = center1, color = 'b')

                Ax[i][j*2].plot(X_axis, Matrix_target2, '-ro', ms = 3)
                Ax[i][j*2].axhline(y = center2, color = 'r')

                Ax[i][j*2].axhline(y = Avg_val, color = 'g')

                if isinstance(Model, np.ndarray): Ax[i][j*2].axhline(y = Model[ncol*i + j], color = 'k')
            
                # Plot horizontal histogram next to optics graph to visualize standard deviation
                hist_data,_,_ = Ax[i][j*2 + 1].hist(Matrix_target1, orientation='horizontal', density=True, color='c', bins=20)
                Ax[i][j*2 + 1].get_yaxis().set_visible(False)
    
    return Center1, Stddev1, Center2, Stddev2

# Determines and Plots Betatron function for IPs
def Betatron_function(s1, s2, sIP, b1, b1_m, b2, b2_m, a1, a1_m, b1_err, b2_err, a1_err, Name1, Name2, title, axis):
    print(title + ':')
    if axis == 'x': axis1 = 'Horizontal'
    else: axis1 = 'Vertical'
    print("%s beta at %s: %f +- %f"%(axis1, Name1, b1, b1_err))
    print("%s alpha at %s: %f +- %f"%(axis1, Name1, a1, a1_err))
    print("%s beta at %s: %f +- %f"%(axis1, Name2, b2, b2_err))
     
    s_star = s1 + a1*b1/(a1**2 + 1)
    beta_star = b1/(a1**2 + 1)

    s_star_m = s1 + a1_m*b1_m/(a1_m**2 + 1)
    beta_star_m = b1_m/(a1_m**2 + 1)

    # Guillaume's Method of calculation s*/b*
    # doesn't use alpha so will be different experimentally
    def fit_beta(x, b1, s1):
        if b1 == 0: return b1
        return b1 + ((x-s1)**2)/b1
    tB = [s1 - sIP, s2 - sIP]
    if np.isnan(b1) or np.isnan(b2): xB = [0, 0]
    else: xB = [b1, b2]
    
    pB = [beta_star, 0]
    parB, covB = curve_fit(fit_beta, tB, xB, p0=pB)
    print("LR values curve fit", parB)

    s_points = [s1, s_star, s2]
    b_points = [b1, beta_star, b2]

    s_points_m = [s1, s_star_m, s2]
    b_points_m = [b1_m, beta_star_m, b2_m]
    
    dbsdb27 = 1/(a1**2 + 1)
    dbsda27 = -2*a1*b1/(a1**2 + 1)**2
    dssdb27 = a1/(a1**2 + 1)
    dssda27 = -(a1**2 - 1)*b1/(a1**2 + 1)**2
    
    b_star_err_prop = np.sqrt((dbsdb27*b1_err)**2 + (dbsda27*a1_err)**2)
    s_star_err_prop = np.sqrt((dssdb27*b1_err)**2 + (dssda27*a1_err)**2)
    
    print("beta star %s: %f +- %f"%(axis, beta_star, b_star_err_prop))
    print("s star %s: %f +- %f"%(axis, s_star, s_star_err_prop))
    print("change in s star %s: %f +- %f"%(axis, s_star - sIP, s_star_err_prop))
    
    err_xpoints = [0, s_star_err_prop, 0]
    err_ypoints = [b1_err, b_star_err_prop, b2_err]

    print()
    print("Madx Values:")
    print("%s beta at %s: %f"%(axis1, Name1, b1_m))
    print("%s alpha at %s: %f"%(axis1, Name1, a1_m))
    print("%s beta at %s: %f"%(axis1, Name2, b2_m))

    tB = [s1 - sIP, s2 - sIP]
    xB = [b1_m, b2_m]
    pB = [beta_star_m, 0]
    parB, covB = curve_fit(fit_beta, tB, xB, p0=pB)
    print("MADx curve fit", parB)

    print("beta star %s: %f"%(axis, beta_star_m))
    print("s star %s: %f"%(axis, s_star_m))
    print("change in s star %s: %f"%(axis, s_star_m - sIP))
    
    # Calculate parabola function from three points
    def calc_parabola_vertex(x1, y1, x2, y2, x3, y3):
        if np.isnan(y1) or np.isnan(y2) or np.isnan(y3): return 0, 0, 0
        denom = (x1-x2) * (x1-x3) * (x2-x3)
        if denom == 0: return 0, 0, 0
        A     = (x3 * (y2-y1) + x2 * (y1-y3) + x1 * (y3-y2)) / denom
        B     = (x3*x3 * (y1-y2) + x2*x2 * (y3-y1) + x1*x1 * (y2-y3)) / denom
        C     = (x2 * x3 * (x2-x3) * y1+x3 * x1 * (x3-x1) * y2+x1 * x2 * (x1-x2) * y3) / denom
        return A, B, C
    
    A, B, C = calc_parabola_vertex(s1, b1, s_star, beta_star, s2, b2)
    A_shift, B_shift, C_shift = calc_parabola_vertex(s1 - s_star, b1, 0, beta_star, s2 - s_star, b2)
    A_m, B_m, C_m = calc_parabola_vertex(s1, b1_m, s_star_m, beta_star_m, s2, b2_m)
    A_m_shift, B_m_shift, C_m_shift = calc_parabola_vertex(s1 - s_star_m, b1_m, 0, beta_star_m, s2 - s_star_m, b2_m)
    
    print()
    print("Experimental Betatron Function:")
    print("%f%s^2 + %d%s + %f"%(A_shift, axis, B_shift, axis, C_shift))

    print("MADx Betatron Function:")
    print("%f%s^2 + %d%s + %f"%(A_m_shift, axis, B_m_shift, axis, C_m_shift))
    print()
    print()

    s_space = np.linspace(s1, s2, 1000)
    Betatron = A*s_space**2 + B*s_space + C
    Betatron_m = A_m*s_space**2 + B_m*s_space + C_m
    
    return s_points, b_points, err_xpoints, err_ypoints, s_points_m, b_points_m, Betatron, Betatron_m

def IP_calculation(Bet_mean, Bet_stdv, Alf_mean, Alf_stdv, Drift_bpm, S_bpms, BetBPMModel, AlfBPMModel, Name_bpms, axis):
    Data = []
    # Set up variables
    def onlyIP(M):
        res = np.zeros(8)
        res[::2] = M[2::6]
        res[1::2] = M[3::6]
        return res
        
    Bet_IP_mean = onlyIP(Bet_mean)
    Bet_IP_stdv = onlyIP(Bet_stdv)
    Alf_IP_mean = onlyIP(Alf_mean)
    Alf_IP_stdv = onlyIP(Alf_stdv)

    IPs = Drift_bpm[1::3]

    sIPs = [0, 639.445027949618, 1277.94839381844, 1917.39342276804]

    # Plot Beta functions at IP 6, 8, 10, 12
    nrow = 2
    ncol = 2
    fig, Ax = plt.subplots(nrow, ncol, figsize=(12, 8), layout='constrained', dpi = 75)
    for i in range(nrow):
        for j in range(ncol):
            index = i*ncol + j
            sIP = sIPs[index]
            C = 3833.845181
            
            # Horizontal calculation
            if axis == "x": title = "Horizontal Betatron function at IP%s"%(6 + index*2)
            else: title = "Vertical Betatron function at IP%s"%(6 + index*2)

            x1, x2 = IPs[index]
            if index == 0: s1 = S_bpms[0][x1] - C
            else: s1 = S_bpms[0][x1]
            s2 = S_bpms[0][x2]
            b1, bx2 = Bet_IP_mean[index*2: index*2 + 2]
            b1_m = BetBPMModel[x1]
            b2_m = BetBPMModel[x2]
            a1 = Alf_IP_mean[index*2]
            a1_m = AlfBPMModel[x1]
            b1_err, b2_err = Bet_IP_stdv[index*2: index*2 + 2]
            a1_err = Alf_IP_stdv[index*2]
            Name1 = Name_bpms[x1]
            Name2 = Name_bpms[x2]
            s_points, b_points, err_xpoints, err_ypoints, s_points_m, b_points_m, Betatron, Betatron_m = Betatron_function(s1, s2, sIP, b1, b1_m, bx2, b2_m, a1, a1_m, b1_err, b2_err, a1_err, Name1, Name2, title, axis)
            Data.append([s_points, err_xpoints, b_points, err_ypoints])

            # Plot Horizontal Beta function
            s_space = np.linspace(s1, s2, 1000)
            Ax[i][j].plot(s_points, b_points, 'o')
            Ax[i][j].plot(s_points_m, b_points_m, 'ok')
            for x, y, text in zip(s_points_m, b_points_m, [Name1, r"$\beta^*_%s$"%axis, Name2]):
                Ax[i][j].annotate(text, # this is the text
                            (x,y), # these are the coordinates to position the label
                            textcoords="offset points", # how to position the text
                            xytext=(0,10), # distance from text to points (x,y)
                            ha='center') # horizontal alignment can be left, right or center
            
            Ax[i][j].plot(s_space, Betatron_m, 'k', label = "MADx")
            Ax[i][j].plot(s_space, Betatron, label = "Experiment")
            Ax[i][j].errorbar(s_points, b_points, xerr = err_xpoints, yerr = err_ypoints, fmt = 'x', ecolor = 'r', capsize = 2)
            Ax[i][j].set_title(title)
            Ax[i][j].set_xlabel("s")
            Ax[i][j].set_ylabel(r"$\beta_%s$(s)"%axis)
            Ax[i][j].legend()

    plt.show()
    return Data

# Main function
def main(args):
    # Parse arguments
    args = argparser(args)
    sdds_file = args.sdds_file
    twiss_file = args.twiss_file
    ring = args.ring
    isNLPlot = args.isNLPlot

    # Initial Data collection
    BPMx, BPMy, Name_xbpms, Name_ybpms, S_xbpms, S_ybpms = data_from_sdds(sdds_file, ring)
    BetxBPMModel, BetyBPMModel, AlfxBPMModel, AlfyBPMModel, BetxModel, AlfxModel, BetyModel, AlfyModel = data_from_twiss(twiss_file)
    start_turn_x, start_turn_y, start_interval, end_interval = initial_plots(BPMx, BPMy)
    tunex, tuney, Bad_xBpms, Bad_yBpms = Tune_Map(BPMx, BPMy)

    # LR Data Analysis
    #  Drift spaces
    if ring == "Blue":
        Drift_bpmx0 = [(163, 165), (167, 0), (2, 4), (23, 25), (27, 28), (30, 32), (53, 55), (57, 60), (62, 64), (83, 84), (86, 87), (89, 90)]
        Drift_bpmy0 = [(162, 164), (166, 0), (2, 4), (24, 26), (28, 29), (31, 33), (53, 54), (56, 59), (61, 63), (83, 84), (86, 87), (89, 90)]
    elif ring == "Yellow":
        Drift_bpmx0 = [(162, 164), (166, 0), (2, 4), 
                      (24, 26), (28, 29), (31, 33), 
                      (53, 55), (57, 60), (62, 63), 
                      (83, 84), (86, 87), (89, 90)]
        Drift_bpmy0 = [(162, 164), (166, 0), (2, 4), 
                      (23, 25), (27, 28), (30, 32), 
                      (53, 55), (57, 60), (62, 63), 
                      (82, 83), (85, 86), (88, 89)]
    else: 
        print("invalid ring (Blue or Yellow)")
        sys.exit()

    # Check if valid drift spaces
    Drift_bpmx = [(None if x in Bad_xBpms else x, None if y in Bad_xBpms else y) for x, y in Drift_bpmx0]
    Drift_bpmy = [(None if x in Bad_yBpms else x, None if y in Bad_yBpms else y) for x, y in Drift_bpmy0]

    # Prepare intervals to look at for LR optics
    IntervalTurns = np.array([i*10 for i in range(start_interval, end_interval + 1)])
    N_optic_turns = len(IntervalTurns)
    IntervalTurnsx = np.array([start_turn_x*np.ones(N_optic_turns), start_turn_x + IntervalTurns], dtype = int).T
    IntervalTurnsy = np.array([start_turn_y*np.ones(N_optic_turns), start_turn_y + IntervalTurns], dtype = int).T

    interval_init = 100
    InitialTurns = np.array([3*i for i in range(N_optic_turns)])
    InitialTurnsx = np.array([start_turn_x + InitialTurns, start_turn_x + InitialTurns + interval_init], dtype = int).T
    InitialTurnsy = np.array([start_turn_y + InitialTurns, start_turn_y + InitialTurns + interval_init], dtype = int).T

    print("\nVary turn interval while keeping starting turn constant: ")
    print("Compute Optics for x:")
    BPMxDriftNames1, Beatx_Matrix1, Betx_Matrix1, Alfx_Matrix1 = Compute_Optics(IntervalTurnsx, Drift_bpmx, Name_xbpms, BetxBPMModel, BPMx, S_xbpms, tunex, ring)

    print("\nCompute Optics for y:")
    BPMyDriftNames1, Beaty_Matrix1, Bety_Matrix1, Alfy_Matrix1 = Compute_Optics(IntervalTurnsy, Drift_bpmy, Name_ybpms, BetyBPMModel, BPMy, S_ybpms, tuney, ring)

    print("\nVary starting turn number while keeping interval constant: ")
    print("Compute Optics for x:")
    BPMxDriftNames2, Beatx_Matrix2, Betx_Matrix2, Alfx_Matrix2 = Compute_Optics(InitialTurnsx, Drift_bpmx, Name_xbpms, BetxBPMModel, BPMx, S_xbpms, tunex, ring)

    print("\nCompute Optics for y:")
    BPMyDriftNames2, Beaty_Matrix2, Bety_Matrix2, Alfy_Matrix2 = Compute_Optics(InitialTurnsy, Drift_bpmy, Name_ybpms, BetyBPMModel, BPMy, S_ybpms, tuney, ring)

    def driftModel(Model, BPMNames):
        ModelDrift = np.zeros(len(BPMNames))
        for i, bpmName in enumerate(BPMNames):
            elem = (bpmName.upper()).replace('-', '_')
            ModelDrift[i] = Model[elem]
        return ModelDrift
    
    # Plot Optics as a function of # Turns used in LR
    print("\nNonlinearity Graphs:\n")
    print("Blue: Vary number of turns; keep initial turn constant")
    print("Red: Vary intital turn; keep number of turns constant")
    print("Green: Average value of their averages")
    print("Black: Model Value")

    name = "Horizontal Beta Beat vs # Turns used in LR"
    Beatx_mean1, Beatx_stdv1, Beatx_mean2, Beatx_stdv2 = NonlinearityGraph(Beatx_Matrix1.T, Beatx_Matrix2.T, BPMxDriftNames1, name, isNLPlot, np.zeros(len(Drift_bpmx)*2))

    name = "Vertical Beta Beat vs # Turns used in LR"
    Beaty_mean1, Beaty_stdv1, Beaty_mean2, Beaty_stdv2 = NonlinearityGraph(Beaty_Matrix1.T, Beaty_Matrix2.T, BPMyDriftNames1, name, isNLPlot, np.zeros(len(Drift_bpmy)*2))

    name = "Horizontal Beta vs # Turns used in LR"
    Betx_mean1, Betx_stdv1, Betx_mean2, Betx_stdv2 = NonlinearityGraph(Betx_Matrix1.T, Betx_Matrix2.T, BPMxDriftNames1, name, isNLPlot, driftModel(BetxModel, BPMxDriftNames1))

    name = "Horizontal Alpha vs # Turns used in LR"
    Alfx_mean1, Alfx_stdv1, Alfx_mean2, Alfx_stdv2 = NonlinearityGraph(Alfx_Matrix1.T, Alfx_Matrix2.T, BPMxDriftNames1, name, isNLPlot, driftModel(AlfxModel, BPMxDriftNames1))

    name = "Vertical Beta vs # Turns used in LR"
    Bety_mean1, Bety_stdv1, Bety_mean2, Bety_stdv2 = NonlinearityGraph(Bety_Matrix1.T, Bety_Matrix2.T, BPMyDriftNames1, name, isNLPlot, driftModel(BetyModel, BPMyDriftNames1))

    name = "Vertical Alpha vs # Turns used in LR"
    Alfy_mean1, Alfy_stdv1, Alfy_mean2, Alfy_stdv2 = NonlinearityGraph(Alfy_Matrix1.T, Alfy_Matrix2.T, BPMyDriftNames1, name, isNLPlot, driftModel(AlfyModel, BPMyDriftNames1))
    if isNLPlot: plt.show()

    # Plot beta functions at IP
    IPDatax1 = IP_calculation(Betx_mean1, Betx_stdv1, Alfx_mean1, Alfx_stdv1, Drift_bpmx, S_xbpms, BetxBPMModel, AlfxBPMModel, Name_xbpms, 'x')

    IPDatay1 = IP_calculation(Bety_mean1, Bety_stdv1, Alfy_mean1, Alfy_stdv1, Drift_bpmy, S_ybpms, BetyBPMModel, AlfyBPMModel, Name_ybpms, 'y')

    IPDatax2 = IP_calculation(Betx_mean2, Betx_stdv2, Alfx_mean2, Alfx_stdv2, Drift_bpmx, S_xbpms, BetxBPMModel, AlfxBPMModel, Name_xbpms, 'x')

    IPDatay2 = IP_calculation(Bety_mean2, Bety_stdv2, Alfy_mean2, Alfy_stdv1, Drift_bpmy, S_ybpms, BetyBPMModel, AlfyBPMModel, Name_ybpms, 'y')

    # Write information to file:
    Datum = IPDatax1, IPDatay1, IPDatax2, IPDatay2
    File_names = ["interval_x", "interval_y", "initial_x", "initial_y"]
    path = os.getcwd()
    dir_path = path + "/IP_results/APEX_7-10-24/"
    isExist = os.path.exists(dir_path)
    if not isExist:
        os.makedirs(dir_path)

    for i in range(len(Datum)):
        file = dir_path + sdds_file + "_results_vary_" + File_names[i]
        if os.path.isfile(file):
            os.remove(file)
        DataFile = open(file,'a',newline='')
        DataWriter = csv.writer(DataFile, dialect='excel')
        for IP_data in Datum[i]:
            for line in IP_data:
                DataWriter.writerow(line)
        DataFile.close()
    
    return 0

if __name__ == "__main__":
    main(sys.argv[1:])