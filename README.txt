README

This script calculates the phase transfer matrix to calculate the linear optics at drift spaces. For each drift space, two bpms were used to calculate the angle coordinates, and a linear regression (LR) was done to find the phase transfer matrix. This method is a data driven method independent of any model; the model is only used in the calculation of beta beat. The locations IR 6, 8, 10, 12 and Magnet sections 3/4 downstream and upstream of IRs are considered. The optics at the IRs are then used to calcuate the betatron function as well as beta* and s*. If one of the bpms is bad at a region, both bpms at that region will not have values since the phase transfer matrix depends on both bpms.

Requirements:
- Need python 3 or more

How to run:
- python3 LR_script.py [sdds_file] [twiss_file] [ring] [isNLPlot]

Arguments:
- sdds_file (required, str): can either be a path to an sdds file or a path to a converted sdds file
- twiss_file (required, str): path to a twiss file
- ring (required, str): Blue or Yellow
- isNLPlot (optional): Plots optics as a function of turn number. Use --no-NLPlot to skip these plots, otherwise leave blank

- Ex1: python3 LR_script.py /operations/app_store/RunData/run_fy23/33916/RHIC/Orbit/TBT/Blue/tbt.Wed_Jun_28_20:13:14_2023.sdds /operations/app_store/Ramps/MADX/Au23-100GeV-e0/Blue/twiss.out_Au23-100GeV-e0 Blue 
- Ex2: python3 LR_script.py data_Wed_Jun_28_20-13-14_2023 twiss.out_Au23-100GeV-e0store Blue --no-NLPlot

Inputs:
- Once the initial set of graphs (First horizontal and vertical BPM and Tunes) are displayed, the user can specify (press enter for default values):
- Starting horizontal turn number for LR (Default 35)
- Starting vertical turn number for LR (Default 535)
- Number of turns to start for Linear Regression (default 100)
- Number of turns to end for Linear Regression (default 400)
- The initial plots can be used as a reference to guide the user what to input

Outputs:
- Mean and standard deviation of first horizontal and vertical BPM
- Initial set of graphs containing tbt data of the first horizontal and vertical BPM
- Tunes horizontal and vertical BPMs
- (Optional:) Plots of the optics as a function of turn number
- Plots of the betatron function according to LR optics and MADx optics
- For every Horizontal/Vertical IP:
	- The average and stdev LR beta value of the bpm downstream and upstream of IP
	- The average and stdev LR alpha value of the bpm downstream of IP
	- LR beta* and s* calculated by curve fit (unused in code, only used to compare to analytical calculation below)
	- LR beta* and s* calculated analytically; stdev calculated using error propagation
	- change in s* is the difference between the LR s* value and the position of the interaction point
	
	- The model beta value of the bpm downstream and upstream of IP
	- The model alpha value of the bpm downstream of IP
	- Model beta* and s* calculated by curve fit (unused in code, only used to compare to analytical calculation below)
	- Model beta* and s* calculated analytically
	- change in s* is the difference between the model s* value and the position of the interaction point
	
	- Experimental formula for the beta function used for plotting
	- Model formula for the beta function used for plotting