import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os
import csv
import LR_script
import subprocess

# Parse Arguments from terminal (README)
# ex: >> python3 IP_results.py twiss.out_pp24-100GeV-e1store Blue --sdds_files 'data_Wed_Jul_10_15_18_01_2024','data_Wed_Jul_10_15_22_53_2024','data_Wed_Jul_10_15_25_35_2024' '--no_g_plot' --model_indices 0,0,1
def argparser(args):
    parser = argparse.ArgumentParser(description='Plot IP results from LR_script.py vs tbt data')
    parser.add_argument('twiss_file', type=str,
                        help='path to RHIC twiss file for optics comparison; ex: /operations/app_store/Ramps/MADx/[beam]/[twiss file]')
    parser.add_argument('ring', type=str,
                        help='Color of Ring used (Blue or Yellow)')
    parser.add_argument('--sdds_files', type=list_of_strings,
                        help='paths to RHIC sdds files (ex: /operations/app_store/RunData/[run year]/[run number]/RHIC/Orbit/TBT/[ring]/[sdds file]) or name of converted sdds file (name of txt file in your directory)')
    parser.add_argument('--model_indices', type=list_of_ints,
                        help='indices of model IP results (ex: --model_indices 0,0,1)')
    parser.add_argument('--no_g_plot', dest='is_g_plot',
                        action='store_false',
                        help="Disable Guillaume's Results; Default is to enable")


    args = parser.parse_args()
    return args

def list_of_strings(arg):
    return arg.split(',')

def list_of_ints(arg):
    return list(map(int, arg.split(',')))

def plot_IP_optics(DataResults, G_Results, Model_Results, Sdds_files, IPs, Name_bpms, axis, is_g_plot):

    # Rearrange the data for plotting:
    s_star_data = DataResults[:, ::4, 1]
    b_star_data = DataResults[:, 2::4, 1]
    b1_data = DataResults[:, 2::4, 0]
    b2_data = DataResults[:, 2::4, 2]

    s_star_data_vary_interval = s_star_data[::2]
    s_star_data_vary_initial = s_star_data[1::2]

    b_star_data_vary_interval = b_star_data[::2]
    b_star_data_vary_initial = b_star_data[1::2]

    b1_data_vary_interval = b1_data[::2]
    b1_data_vary_initial = b1_data[1::2]

    b2_data_vary_interval = b2_data[::2]
    b2_data_vary_initial = b2_data[1::2]

    # Plot results over runs
    nrow = 4 
    ncol = 4

    fig, Ax = plt.subplots(nrow, ncol, figsize=(15, 8), layout='constrained')
    if axis == "x": fig.suptitle("Horizontal Beta from different TBT Data of APEX Experiment")
    else: fig.suptitle("Vertical Beta from different TBT Data of APEX Experiment")

    x_axis = np.linspace(1, len(Sdds_files), len(Sdds_files))
    for i, file in enumerate(Sdds_files):
        print("%d:"%(i + 1), file)
    print()

    print("blue: vary number of turns in an interval")
    print("red: vary initial turn number")
    print("green: average value of blue and red")
    sIPs = [0, 639.445027949618, 1277.94839381844, 1917.39342276804]
    for i in range(nrow):
        x1, x2 = IPs[i]
        Name1 = Name_bpms[x1]
        Name2 = Name_bpms[x2]
        
        Ax[i][0].plot(x_axis, b1_data_vary_interval[:, i], '-ob', ms = 3)
        Ax[i][1].plot(x_axis, b_star_data_vary_interval[:, i], '-ob', ms = 3)
        Ax[i][2].plot(x_axis, s_star_data_vary_interval[:, i] - sIPs[i], '-ob', ms = 3)
        Ax[i][3].plot(x_axis, b2_data_vary_interval[:, i], '-ob', ms = 3)

        Ax[i][0].plot(x_axis, b1_data_vary_initial[:, i], '-or', ms = 3)
        Ax[i][1].plot(x_axis, b_star_data_vary_initial[:, i], '-or', ms = 3)
        Ax[i][2].plot(x_axis, s_star_data_vary_initial[:, i] - sIPs[i], '-or', ms = 3)
        Ax[i][3].plot(x_axis, b2_data_vary_initial[:, i], '-or', ms = 3)

        Ax[i][0].plot(x_axis, 
                      (b1_data_vary_interval[:, i] + b1_data_vary_initial[:, i])/2, '-og', ms = 3)
        Ax[i][1].plot(x_axis, (b_star_data_vary_interval[:, i] + b_star_data_vary_initial[:, i])/2, '-og', ms = 3)
        Ax[i][2].plot(x_axis, (s_star_data_vary_interval[:, i] + s_star_data_vary_initial[:, i])/2 - sIPs[i], '-og', ms = 3)
        Ax[i][3].plot(x_axis, (b2_data_vary_interval[:, i] + b2_data_vary_initial[:, i])/2, '-og', ms = 3)

        if is_g_plot:
            if i == 0: print("magenta: Guillaume's Results")
            Ax[i][0].plot(x_axis, G_Results.T[i, 0], '-om', ms = 3)
            Ax[i][1].plot(x_axis, G_Results.T[i, 1], '-om', ms = 3)
            Ax[i][2].plot(x_axis, G_Results.T[i, 2], '-om', ms = 3)
            Ax[i][3].plot(x_axis, G_Results.T[i, 3], '-om', ms = 3)

        if i == 0: print("black: Results from ORM data")
        Ax[i][0].plot(x_axis, Model_Results[i][:, 0], '-ok', ms = 3)
        Ax[i][1].plot(x_axis, Model_Results[i][:, 1], '-ok', ms = 3)
        Ax[i][2].plot(x_axis, Model_Results[i][:, 2], '-ok', ms = 3)
        Ax[i][3].plot(x_axis, Model_Results[i][:, 3], '-ok', ms = 3)

        Ax[i][0].set_title(Name1) 
        Ax[i][0].set_ylabel(r"$\beta [m]$") 
        
        Ax[i][1].set_title(r"$\beta^*_%s$"%axis) 
        Ax[i][1].set_ylabel(r"$\beta$ [m]") 
        
        Ax[i][2].set_title(r"$\Delta s^*_%s$"%axis)
        Ax[i][2].set_ylabel("s [m]") 
        
        Ax[i][3].set_title(Name2)
        Ax[i][3].set_ylabel(r"$\beta$ [m]") 

    print()
    print()

    plt.show()

def g_Data(Filenames, Name_xbpms, Name_ybpms, S_xbpms, S_ybpms):
    Gx_Data = []
    Gy_Data = []
    Namex_ring = []
    Namey_ring = []
    sx_ring = []
    sy_ring = []
    Betx_ring = []
    Bety_ring = []
    for i, filename in enumerate(Filenames):
        IP_vals = []
        print(filename)
        with open(filename) as output:
            Lines = csv.reader(output, delimiter = ',', quotechar = '|')
            Lines_list = list(Lines)
            Namex_list = Lines_list[0]
            sx_list = np.array(list(map(float, Lines_list[1])))
            Betx_list = np.array(list(map(float, Lines_list[2])))**2
            Namey_list = Lines_list[3]
            sy_list = np.array(list(map(float, Lines_list[4])))
            Bety_list = np.array(list(map(float, Lines_list[5])))**2
            for line_index, line in enumerate(Lines_list[6:-4]):
                IP_vals.append(np.array(list(map(float, line))))

        IP_vals = np.array(IP_vals)

        Namex_ring.append(Namex_list)
        Namey_ring.append(Namey_list)
        sx_ring.append(sx_list)
        sy_ring.append(sy_list)
        Betx_ring.append(Betx_list)
        Bety_ring.append(Bety_list)
        
        IP_search = ['g5-bx', 'g6-bx', 'g7-bx', 'g8-bx', 'g9-bx', 'g10-bx', 'g11-bx', 'g12-bx']
        IPx_indices = np.zeros(len(IP_search), dtype = int)
        IPy_indices = np.zeros(len(IP_search), dtype = int)
        for i, ip in enumerate(IP_search):
            IPx_indices[i] = Namex_list.index(ip)
            IPy_indices[i] = Namey_list.index(ip)
            
        Betx_IP_list = Betx_list[IPx_indices]
        sx_IP_list = sx_list[IPx_indices]
        C = 3833.845181
        sx_IP_list[0] -= C
        
        Bety_IP_list = Bety_list[IPy_indices]
        sy_IP_list = sy_list[IPy_indices]
        sy_IP_list[0] -= C

        gx_data = np.array([Betx_IP_list[::2],
                    IP_vals[::2, 0],
                    IP_vals[::2, 1],
                    Betx_IP_list[1::2]])
        gy_data = np.array([Bety_IP_list[::2],
                    IP_vals[1::2, 0],
                    IP_vals[1::2, 1],
                    Bety_IP_list[1::2]])
        
        Gx_Data.append(gx_data)
        Gy_Data.append(gy_data)

    # Plot Beta around Ring
    def Mod_Bet_ring(Bet_ring, Name_ring, Name_bpms):
        Mod_Bet_ring = []
        for i in range(len(Bet_ring)):
            mod_bet_ring = np.zeros(len(Name_bpms))
            j = 0
            for bpm in Name_bpms:
                if bpm in Name_ring[i]:
                    # print(Betx_ring[i], Name_ring[i].index(bpm))
                    mod_bet_ring[j] = Bet_ring[i][Name_ring[i].index(bpm)]
                j += 1
            Mod_Bet_ring.append(mod_bet_ring)
        return np.array(Mod_Bet_ring)

    Mod_Betx_ring = Mod_Bet_ring(Betx_ring, Namex_ring, Name_xbpms)
    Mod_Bety_ring = Mod_Bet_ring(Bety_ring, Namey_ring, Name_ybpms)
    
    nrow = 2
    ncol = 1
    fig, Ax = plt.subplots(nrow, ncol, figsize=(14, 8), layout='constrained')
    fig.suptitle("Beta around Ring")
    Ax[0].plot(S_xbpms[0], Mod_Betx_ring[0], "-k", label = "First")
    if len(Mod_Betx_ring) > 1: Ax[0].plot(S_xbpms[0], Mod_Betx_ring[-2], label = "Prev", alpha = .75)
    Ax[0].plot(S_xbpms[0], Mod_Betx_ring[-1], label = "Curr", alpha = .5)
    Ax[1].set_ylabel(r"$\beta_x [m]$")
    Ax[0].legend()

    Ax[1].plot(S_ybpms[0], Mod_Bety_ring[0], "-k", label = "First")
    if len(Mod_Bety_ring) > 1: Ax[1].plot(S_ybpms[0], Mod_Bety_ring[-2], label = "Prev", alpha = .75)
    Ax[1].plot(S_ybpms[0], Mod_Bety_ring[-1], label = "Curr", alpha = .5)
    Ax[1].set_ylabel(r"$\beta_y [m]$")
    Ax[1].set_xlabel("s")
    Ax[1].legend()

    plt.show()

    return np.array(Gx_Data), np.array(Gy_Data)

def get_values(file):
    Values = []
    with open("../sstar-wire/files_e1/%s"%file) as output:
        Lines = csv.reader(output, delimiter = ',', quotechar = '|')
        for i, row in enumerate(Lines):
            Values.append(list(map(float, row)))
    return np.array(Values)

def get_model_data(Model_indices):
    APEX_dir = "/APEX_7-10-24"
    # APEX_dir = "/APEX_8-21-24_pos"
    # APEX_dir = "/APEX_8-21-24_neg"
    IP_data = []
    for i in range(4):
        IP_data_file = APEX_dir + "/IP%d_data.csv"%(2*i + 6)
        IP_data.append(get_values(IP_data_file))

    IPx_data_plot = []
    IPy_data_plot = []
    for i in range(4):
        IPx_split = []
        IPy_split = []
        for j in Model_indices:
            IPx_split.append(IP_data[i][j][:4])
            IPy_split.append(IP_data[i][j][4:])
        IPx_data_plot.append(np.array(IPx_split))
        IPy_data_plot.append(np.array(IPy_split))

    return np.array(IPx_data_plot), np.array(IPy_data_plot)

def main(args):
    argparser(args)
    args = argparser(args)
    Sdds_files = args.sdds_files
    Model_indices = args.model_indices
    twiss_file = args.twiss_file
    ring = args.ring
    is_g_plot = args.is_g_plot

    File_names = ["interval_x", "interval_y", "initial_x", 
                  "initial_y"]
    # Sdds_files = ['data_Wed_Jul_10_15_18_01_2024',
    #               'data_Wed_Jul_10_15_22_53_2024',
    #               'data_Wed_Jul_10_15_25_35_2024',]
    #             #   'data_Wed_Jul_10_15_34_02_2024']
    # twiss_file = 'twiss.out_pp24-100GeV-e1store'
    # ring = "Blue"

    path = os.getcwd()
    DataResultsx = []
    DataResultsy = []
    g_Filenames = []
    for i, sdds_file in enumerate(Sdds_files):
        # Convert sdds files if needed
        if sdds_file[-5:] == ".sdds":
            index_ring = sdds_file.find(ring)
            index_sdds = sdds_file.find('.sdds')
            sdds_file = "data_" + sdds_file[index_ring: index_sdds].replace(":", "_")
            if not os.path.exists(sdds_file):
                subprocess.run(["sddsconvert", "-ascii", sdds_file, sdds_file])
        # Check if IP results are already produced for the given files
        if not os.path.exists(str(path) + "/IP_results/" + sdds_file + "_results_vary_" + File_names[0]): 
            LR_script.main([sdds_file, twiss_file, ring, "--no-NLPlot"])
        g_Filenames.append("g_IP_results/" + sdds_file[5:])
    
        # Load Data from IP Results
        def Load_data(file):
            Data = []
            with open(file) as output:
                Lines = csv.reader(output, delimiter = ',', quotechar = '|')
                for line in Lines:
                    Data.append(list(map(float, line)))
            return Data
        for j, file_name in enumerate(File_names):
            file = "IP_results/" + sdds_file + "_results_vary_" + file_name
            if j%2 == 0: DataResultsx.append(Load_data(file))
            else: DataResultsy.append(Load_data(file))
    DataResultsx = np.array(DataResultsx)
    DataResultsy = np.array(DataResultsy)

    # Blue         
    xIPs = [(167, 0), (27, 28), (57, 60), (86, 87)]
    yIPs = [(166, 0), (28, 29), (56, 59), (86, 87)]
    
    BPMx, BPMy, Name_xbpms, Name_ybpms, S_xbpms, S_ybpms = LR_script.data_from_sdds(Sdds_files[0], ring)

    IPx_model_res, IPy_model_res = get_model_data(Model_indices)

    Gx_Data, Gy_Data = g_Data(g_Filenames, Name_xbpms, Name_ybpms, S_xbpms, S_ybpms)
    plot_IP_optics(DataResultsx, Gx_Data, IPx_model_res, Sdds_files, xIPs, Name_xbpms, "x", is_g_plot)
    plot_IP_optics(DataResultsy, Gy_Data, IPy_model_res, Sdds_files, yIPs, Name_ybpms, "y", is_g_plot)

    return 0

if __name__ == "__main__":
    main(sys.argv[1:])