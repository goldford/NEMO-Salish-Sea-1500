import scipy.io
import os
import pickle
import pandas as pd
import numpy as np
import scipy as sp
import scipy.stats
import sys
from matplotlib import rcParams
from matplotlib.markers import MarkerStyle
import matplotlib
import matplotlib.ticker as ticker
import matplotlib.colors as clr
from matplotlib.ticker import ScalarFormatter
import math
from math import log10, floor
from array import array
import numbers
from matplotlib.lines import Line2D
import warnings
from matplotlib.markers import MarkerStyle
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#Taylor axes code

from matplotlib.projections import PolarAxes
import mpl_toolkits.axisartist.floating_axes as FA
import mpl_toolkits.axisartist.grid_finder as GF

def get_taylor_diagram_axes(fig, rect, refstd, srange, contour_levs, extend=False):
    
    tr = PolarAxes.PolarTransform()
    
    # Correlation labels
    rlocs = np.array([0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1])
    if extend:
        # Diagram extended to negative correlations
        tmax = np.pi
        rlocs = np.concatenate((-rlocs[:0:-1], rlocs))
    else:
        # Diagram limited to positive correlations
        tmax = np.pi/2
    
    tlocs = np.arccos(rlocs)        # Conversion to polar angles
    gl1 = GF.FixedLocator(tlocs)    # Positions
    tf1 = GF.DictFormatter(dict(zip(tlocs, map(str, rlocs))))

    # Standard deviation axis extent (in units of reference stddev)
    smin = srange[0] * refstd
    smax = srange[1] * refstd
    
    
    ghelper = FA.GridHelperCurveLinear(
        tr,
        extremes=(0, tmax, smin, smax),
        grid_locator1=gl1, 
        tick_formatter1=tf1
    )

    ax = FA.FloatingSubplot(fig, rect, grid_helper=ghelper)
    
    fig.add_subplot(ax)

    # Adjust axes
    ax.axis["top"].set_axis_direction("bottom")   # "Angle axis"
    ax.axis["top"].toggle(ticklabels=True, label=True)
    ax.axis["top"].major_ticklabels.set_axis_direction("top")
    ax.axis["top"].label.set_axis_direction("top")
    ax.axis["top"].label.set_text("Correlation")

    ax.axis["left"].set_axis_direction("bottom")  # "X axis"
    ax.axis["left"].label.set_text("Standard deviation")

    ax.axis["right"].set_axis_direction("top")    # "Y-axis"
    ax.axis["right"].toggle(ticklabels=True)
    ax.axis["right"].major_ticklabels.set_axis_direction(
        "bottom" if extend else "left")

    if smin:
        ax.axis["bottom"].toggle(ticklabels=False, label=False)
    else:
        ax.axis["bottom"].set_visible(False)          # Unused

    # Graphical axes
    _ax = ax
    ax = ax.get_aux_axes(tr)   # Polar coordinates

    # Add reference point and stddev contour
    l, = ax.plot([0], refstd, 'k*', ls='', ms=10, label="label")

    t = np.linspace(0, tmax)
    r = np.zeros_like(t) + refstd
    ax.plot(t, r, 'k--', label='_')
    
    rs, ts = np.meshgrid(np.linspace(smin, smax), np.linspace(0, tmax))
    
    # Compute centered RMS difference
    rms = np.sqrt(refstd**2 + rs**2 - 2*refstd*rs*np.cos(ts))

    
    _ax.grid(ls='-',color='#969696',lw=0.8)
    _ax.axis[:].major_ticks.set_tick_out(True) 
    
    
    _ax.axis[:].major_ticks.set_tick_out(True)
    _ax.xaxis.set_tick_params(labelsize=2)
    _ax.axis["top"].label.set_text("Correlation")
    _ax.axis["left"].label.set_text("Normalized standard deviation")
    
    
    
    contours = ax.contour(ts, rs, rms, colors='0.6', levels = contour_levs, linewidths = 1, linestyles='--')
    plt.clabel(contours, inline=1, fontsize=10, fmt='%.2f', colors='#5A5A5A')

#     # Collect sample points for latter use (e.g. legend)
#     self.samplePoints = [l]   

    return ax


def plot_taylor_data(ax,modelruns_info,var,augment_rmse,buoy_list,o_markers,o_markerlabel):
    # get data and plot
    for run in modelruns_info.keys():

        if var == "temperature":
            print("Taylor SST", run)
        elif var == "salinity":
            print("Taylor SSS", run)
        else:
            print("No variable chosen")
            break

        path1 = modelruns_info[run]['path']
        shortcode = modelruns_info[run]['shortcode']
        colour_m = modelruns_info[run]['colour']

        # plot markers
        facecolor1 = colour_m
        facecolor2 = colour_m
        edgecolor1 = colour_m
        edgecolor2 = colour_m
        medgew1 = 0.2
        medgew2 = 0.2

        # get lighthouse data
        LH_pic_data = 'LH_class4_' + shortcode + '_hindcast.pickle'
        LH_scores = pickle.load(open(os.path.join(path1,LH_pic_data), 'rb'))             

        # to do - am returning data for target unnecessarily - modify get_lh_stats into 2 functions
        lh_stats_taylor_np, lh_stats_target_np = get_lh_stats(LH_scores, run, augment_rmse, var)

        if var == "temperature":
            # SST buoy data
            SSTbuoy_pic_data = "SST_class4_" + run + ".pickle"
            SSTbuoy_scores = pickle.load(open(os.path.join(path1,SSTbuoy_pic_data), 'rb'))

            # looks like:
            # SSTbuoy_scores['SalishSea1500-RUN203']['ptaw1']['filt_interp_data']['scores'].keys()
            # dict_keys(['skill1981', 'mean_obs', 'mean_mod', 'bias', 'crmse', 'gamma2', 'rmse', 'mae', 'mad', 'pearson', 'stdev_obs', 'stdev_mod'])
            buoy_stats_taylor_np, buoy_stats_target_np = get_SSTbuoy_stats(SSTbuoy_scores, run, augment_rmse, buoy_list)
            SST_taylor_np_buoys_and_LH = np.concatenate((lh_stats_taylor_np,buoy_stats_taylor_np))

            ccoef = SST_taylor_np_buoys_and_LH[:,6].astype('f')
            sdev = SST_taylor_np_buoys_and_LH[:,3].astype('f')

        elif var == "salinity":
            lh_stats_taylor_S_np, lh_stats_target_S_np = get_lh_stats(LH_scores, run, augment_rmse, var)

            ccoef = lh_stats_taylor_S_np[:,6].astype('f')
            sdev = lh_stats_taylor_S_np[:,3].astype('f')

        pairs_stats = []
        for i in range(len(sdev+1)):
            if i == 0:
                stn = "reference"
                sig = 1
                cc = 1
            else:
                stn = o_markerlabel[i-1]
                sig = sdev[i-1]
                cc = ccoef[i-1]
            pairs_stats.append([stn,[sig,cc]])

        # Add the model markers to Taylor diagram
        i = 0
        for i, (stn, stats) in enumerate(pairs_stats):
            # avoid douvble-plotting ref
            if stn == 'reference':
                continue
            stddev = stats[0]
            corrcoef = stats[1]
            #print(corrcoef)

            ax.plot(np.arccos(corrcoef), 
                    stddev, 
                    marker=o_markers[stn]['symbol'],
                    ms=o_markers[stn]['size'],
                    ls='',
                    #mfc=colors[i],
                    #mfc=o_markers[stn]['facecolor'],
                    mfc=colour_m,
                    #mec=colors[i],
                    mec=edgecolor1,
                    mew=medgew1,
                    #label="%d" % (i+1)),
                    label=o_markerlabel[i]
                   )
            

    return ax
            
            
def plot_target_data1(ax,modelruns_info,var,augment_rmse, options):


    # loop over model runs, add data
    for run in modelruns_info.keys():
      
        if var == "temperature":
            print("Target SST", run)
        elif var == "salinity":
            print("Target SSS", run)
        else:
            print("No variable chosen")
            break
    
        path1 = modelruns_info[run]['path']
        shortcode = modelruns_info[run]['shortcode']
        colour_m = modelruns_info[run]['colour']
        
        # markers
        facecolor1 = colour_m
        facecolor2 = colour_m
        edgecolor1 = colour_m
        edgecolor2 = colour_m
        medgew1 = 0.2
        medgew2 = 0.2
        

        # get lighthouse data
        LH_pic_data = 'LH_class4_' + shortcode + '_hindcast.pickle'
        LH_scores = pickle.load(open(os.path.join(path1,LH_pic_data), 'rb'))             
        lh_stats_taylor_np, lh_stats_target_np = get_lh_stats(LH_scores, run, augment_rmse, var)
        
        # if temp, then use buoys otherwise just lh data
        if var == "temperature":
             # SST buoy data
            SSTbuoy_pic_data = "SST_class4_" + run + ".pickle"
            #print(os.path.join(path1,SSTbuoy_pic_data)
            SSTbuoy_scores = pickle.load(open(os.path.join(path1,SSTbuoy_pic_data), 'rb'))

            # looks like:
            # SSTbuoy_scores['SalishSea1500-RUN203']['ptaw1']['filt_interp_data']['scores'].keys()
            # dict_keys(['skill1981', 'mean_obs', 'mean_mod', 'bias', 'crmse', 'gamma2', 'rmse', 'mae', 'mad', 'pearson', 'stdev_obs', 'stdev_mod'])
            buoy_stats_taylor_np, buoy_stats_target_np = get_SSTbuoy_stats(SSTbuoy_scores, run, augment_rmse, buoy_list)
            SST_taylor_np_buoys_and_LH = np.concatenate((lh_stats_taylor_np,buoy_stats_taylor_np))
            SST_target_np_buoys_and_LH = np.concatenate((lh_stats_target_np,buoy_stats_target_np))

            Bs = SST_target_np_buoys_and_LH[:,5].astype('f')
            RMSDs =  SST_target_np_buoys_and_LH[:,4].astype('f')

        elif var == "salinity":
            Bs = lh_stats_target_np[:,5].astype('f')
            RMSDs =  lh_stats_target_np[:,4].astype('f')

        ax = plot_target_data2(ax,RMSDs,Bs,options, facecolor1, edgecolor1, medgew1)
        
    return ax

def plot_target_data2(ax,RMSDs,Bs,options, facecolor1, edgecolor1, medgew1):

    x = RMSDs
    y = Bs 
    
    # Set font and marker size
    fontSize = matplotlib.rcParams.get('font.size') 

    o_markerlabel = options['o_markerlabel']
    o_markers = options['o_markers']
    o_markerstyle = options['o_markerstyle'] # added by GO - "font" / "symbol"
    o_markersize = options['o_markersize']
    o_markersymbol = options['o_markersymbol'] #'+','o','x','s','d','^','v','p','h','*'
    o_markercolor = options['o_markercolor'] # single color to use for all markers (Default: None)
    o_markercolors = options['o_markercolors'] # if None, considers 'markercolor' only; dictionary with two colors as keys ('face', 'edge')
#                             or None. If None or 'markerlegend' == 'on' then
#                             considers only the value of 'markercolor'. (Default: None)
    o_markerlabelcolor = options['o_markerlabelcolor']
    markerdisplayed = options['markerdisplayed']
    markerobs = options['markerobs']
    markerlegend = options['markerlegend']
    alpha = options['alpha']
    stylebias = options['stylebias']
    numberpanels = options['numberpanels']
    face = options['face']
    edge = options['edge']
    cmap = options['cmap']
    cmap_vmin = options['cmap_vmin']
    cmap_vmax = options['cmap_vmax']
    cmap_marker = options['cmap_marker']
    cmapzdata = options['cmapzdata']
    colframe = options['colframe']
    colormap = options['colormap']
    labelweight = options['labelweight']
    locationcolorbar = options['locationcolorbar']
    axismax = options['axismax']
    
    target_options_file = options['target_options_file']
    titlecolorbar = options['titlecolorbar']

    # Check enough labels provided if markerlabel provided. Not a problem if labels
    # provided via the markers option.
    numberLabel = len(o_markerlabel)
    if numberLabel > 0:
        if isinstance(o_markerlabel, list) and numberLabel < len(RMSDs):
            raise ValueError('Insufficient number of marker labels provided.\n' +
                             'target: No. labels=' + str(numberLabel) + ' < No. markers=' +
                             str(len(RMSDs)) + '\n' +
                             'taylor: No. labels=' + str(numberLabel+1) + ' < No. markers=' +
                             str(len(RMSDs)+1))
        elif isinstance(o_markerlabel, dict) and numberLabel > 70: # what is this? GLO why 70?
            raise ValueError('Insufficient number of marker labels provided.\n' +
                             'target: No. labels=' + str(numberLabel) + ' > No. markers= 70')

    if markerlegend == 'on':
        # Check that marker labels have been provided
        if o_markerlabel == '' and o_markers == None:
            raise ValueError('No marker labels provided.')

        # Plot markers of different color and symbols with labels displayed in a legend
        limit = axismax
        hp = ()
        rgba = None

        if o_markers is None:

            # Define default markers (function)
    #         marker, markercolor = get_default_markers(X, option)
            #get_default_markers:
    #         alpha = option['alpha']

            # Define list of marker symbols and colros
            kind = ['+','o','x','s','d','^','v','p','h','*']
            colorm = ['r','b','g','c','m','y','k','gray']
            if len(RMSDs) > 80:
                print('You must introduce new markers to plot more than 70 cases.')
                print('The ''marker'' character array need to be extended inside the code.')

            if len(RMSDs) <= len(kind):
                # Define markers with specified color
                marker = []
                markercolor = []
                if o_markercolor is None:
                    for i, color in enumerate(colorm):
                        rgba = clr.to_rgb(color) + (alpha,)
                        marker.append(kind[i] + color)
                        markercolor.append(rgba)
                else:
                    rgba = clr.to_rgb(o_markercolor) + (alpha,)
                    for symbol in kind:
                        marker.append(symbol + o_markercolor)
                        markercolor.append(rgba)
            else:
                # Define markers and colors using predefined list
                marker = []
                markercolor = []
                for color in colorm:
                    for symbol in kind:
                        marker.append(symbol + color)
                        rgba = clr.to_rgb(color) + (alpha,)
                        markercolor.append(rgba)
            # END get_default_markers


            # Plot markers at data points
            labelcolor = []
            markerlabel = []
            for i, xval in enumerate(RMSDs):
                if abs(RMSDs[i]) <= limit and abs(Bs[i]) <= limit:
                    h = ax.plot(RMSDs[i],
                                Bs[i],
                                marker[i], 
                                markersize = o_markersize,
                                markerfacecolor = markercolor[i],
                                markeredgecolor = markercolor[i][0:3] + (1.0,),
                                markeredgewidth = 1)
                    hp += tuple(h)
                    labelcolor.append(o_markerlabelcolor)
                    markerlabel.append(o_markerlabel[i])

        # if there's an o_markers dictionary
        else:
            # Obtain markers from option['markers']
            #labels, labelcolor, marker, markersize, markerfacecolor, markeredgecolor = \
                #get_single_markers(option['markers'])

            #get_single_markers -->
            if o_markers is None:
                raise ValueError("Empty dictionary provided for option['markers']")

            labelcolor = []
            marker = []
            markerfacecolor = []
            markeredgecolor = []
            markerlabel = []
            markersize = []
            markeredgewidth = []

    #         if o_markerstyle == "font":
    #             marker_style.update(markeredgecolor="none", markersize=15)

            # Iterate through keys in dictionary
            for key in o_markers:
#                 color = o_markers[key]['facecolor']
                color = facecolor1
                symbol = o_markers[key]['symbol']

                if o_markerstyle != "font":
                    SymbolColor = symbol + color
                    marker.append(SymbolColor)
                else:
                    marker.append(symbol)
                    fontFamily = rcParams.get('font.family')

                markersize.append(o_markers[key]['size'])
#                 markerfacecolor.append(color)
                markerfacecolor.append(facecolor1)
#                 markeredgecolor.append(o_markers[key]['edgecolor'])
                markeredgecolor.append(edgecolor1)
#                 markeredgewidth.append(o_markers[key]['edgewidth'])
                markeredgewidth.append(medgew1)
                markerlabel.append(key) # store label
#                 labelcolor.append(o_markers[key]['labelcolor'])
                labelcolor.append(facecolor1)
                
            # end get_single_markers

            # Plot markers at data points
            for i, xval in enumerate(RMSDs):
                if abs(RMSDs[i]) <= limit and abs(Bs[i]) <= limit:
                    h = ax.plot(RMSDs[i],
                                Bs[i],
                                marker=marker[i],
                                markersize = markersize[i],
                                markerfacecolor = markerfacecolor[i],
                                markeredgecolor = markeredgecolor[i],
                                markeredgewidth = markeredgewidth[i])
                    hp += tuple(h)
                    #markerlabel.append(labels[i])

    else:
        # Plot markers as dots of a single color with accompanying labels
        limit = axismax

        # Define edge and face colors of the markers
        #                                             default_key, dict_key, key_key
        #edge_color = get_from_dict_or_default(option, 'markercolor', 'markercolors', 'edge')
        # get_from_dict_or_default ->
        if o_markercolors is None:
            edge_color = o_markercolor
        elif edge is None:
            edge_color = o_markercolor
        elif edge != None:
            edge_color = edge
        # end get_from_dict_or_default

        if edge_color is None: edge_color = 'r'

        #face_color = get_from_dict_or_default(option, 'markercolor', 'markercolors', 'face')
        # get_from_dict_or_default ->
        if o_markercolors is None:
            face_color = o_markercolor
        elif face is None:
            face_color = o_markercolor
        elif face != None:
            face_color = face
        # end get_from_dict_or_default

        if face_color is None: face_color = edge_color
        face_color = clr.to_rgb(face_color) + (alpha,)

        labelcolor = []
        for i in range(len(RMSDs)):
            xval, yval = RMSDs[i], Bs[i]
            if abs(xval) <= limit and abs(yval) <= limit:
                ax.plot(xval, yval, o_markersymbol,markersize=o_markersize, markerfacecolor=face_color,markeredgecolor=edge_color)
                labelcolor.append(o_markerlabelcolor)

                # Check if marker labels provided
                if type(o_markerlabel) is list:
                    # Label marker
                    ax.text(xval, yval, o_markerlabel[i],color=o_markerlabelcolor,verticalalignment='bottom',horizontalalignment='right',fontsize=fontSize)

            del i, xval, yval
            
    return ax


# functions to load LH and buoy data
def get_lh_stats(LH_scores, modelrun, augment_rmse, var):
    
    lh_stats_taylor = []
    lh_stats_target = []

    for lh in LH_scores[modelrun].keys():
        
        pears = LH_scores[modelrun][lh]['filt_interp_data']['scores'][var]['pearson']
        mod_stdev = LH_scores[modelrun][lh]['filt_interp_data']['scores'][var]['stdev_mod']
        obs_stdev = LH_scores[modelrun][lh]['filt_interp_data']['scores'][var]['stdev_obs']

        mod_norm_stdev = mod_stdev / obs_stdev

        crmse = LH_scores[modelrun][lh]['filt_interp_data']['scores'][var]['crmse']

        ncrmse = crmse / obs_stdev   

        # stats needed for taylor
        lh_stats_taylor.append([lh,
                                obs_stdev,
                                mod_stdev,
                                mod_norm_stdev,
                                crmse,
                                ncrmse,
                                pears])

        # stats needed for target
        rmse = LH_scores[modelrun][lh]['filt_interp_data']['scores'][var]['rmse']
        bias = LH_scores[modelrun][lh]['filt_interp_data']['scores'][var]['bias']
        nrmse = rmse / obs_stdev

        if augment_rmse == True:
            if (mod_stdev - obs_stdev) < 0:
                rmse = rmse * -1
                nrmse = nrmse * -1
                crmse = crmse * -1
                ncrmse = ncrmse * -1

        lh_stats_target.append([lh,
                                rmse,
                                crmse,
                                nrmse,
                                ncrmse,
                                bias
                               ])

    lh_stats_taylor_np = np.array(lh_stats_taylor)
    lh_stats_target_np = np.array(lh_stats_target)
    
    
    return(lh_stats_taylor_np, lh_stats_target_np)


def get_SSTbuoy_stats (SSTbuoy_scores, modelrun, augment_rmse, buoy_list):
    buoy_stats_taylor = []
    buoy_stats_target = []

    for buoy in SSTbuoy_scores[modelrun].keys():

        if buoy == "C46134_2001-02_2016-12.nc":
            buoy_bettername = '46134'
        elif buoy == "C46131_1992-10_2021-01.nc":
            buoy_bettername = '46131'
        elif buoy == "C46146_1992-03_2022-06.nc":
            buoy_bettername = '46146'
        elif buoy == "C46182_1989-09_1991-11.nc":
            buoy_bettername = '46182'
        else:
            buoy_bettername = buoy
        
        if (buoy_bettername not in buoy_list):
            continue
        
        obs_stdev = SSTbuoy_scores[modelrun][buoy]['filt_interp_data']['scores']['stdev_obs']
        mod_stdev = SSTbuoy_scores[modelrun][buoy]['filt_interp_data']['scores']['stdev_mod']
        mod_norm_stdev = mod_stdev / obs_stdev
        crmse = SSTbuoy_scores[modelrun][buoy]['filt_interp_data']['scores']['crmse']
        ncrmse = crmse / obs_stdev
        pears = SSTbuoy_scores[modelrun][buoy]['filt_interp_data']['scores']['pearson']
        
        # for taylors
        buoy_stats_taylor.append([buoy_bettername,
                                  obs_stdev,
                                  mod_stdev,
                                  mod_norm_stdev,
                                  crmse,
                                  ncrmse,
                                  pears]) 
        # for target
        rmse = SSTbuoy_scores[modelrun][buoy]['filt_interp_data']['scores']['rmse']
        bias = SSTbuoy_scores[modelrun][buoy]['filt_interp_data']['scores']['bias']
        nrmse = rmse / obs_stdev

        if augment_rmse == True:
            if (mod_stdev - obs_stdev) < 0:
                rmse = rmse * -1
                crmse = crmse * -1
                nrmse = nrmse * -1
                ncrmse = ncrmse * -1

        buoy_stats_target.append([buoy_bettername,
                                  rmse,
                                  crmse,
                                  nrmse,
                                  ncrmse,
                                  bias])

    # np is easier
    buoy_stats_taylor_np = np.array(buoy_stats_taylor)
    buoy_stats_target_np = np.array(buoy_stats_target)
    
    return(buoy_stats_taylor_np, buoy_stats_target_np)



### TARGET figures, diagrams
# SM code pasted below (instead of difficult to follow separate .py)

# import sys
# !{sys.executable} -m pip3 install skillmetrics

# functions required for SkillMetrics
# functions
def find_exp(number) -> int:
    base10 = log10(abs(number))
    return floor(base10)

def use_sci_notation(value):
    '''
    Boolean function to determine if scientific notation to be used for value
 
    Input:
    Value : value to be tested
 
    Return:
        True if absolute value is > 100 or < 1.e-3
        False otherwise

    Author: Peter A. Rochford
        Symplectic, LLC

    Created on May 10, 2022
    '''
    if (abs(value)>0 and abs(value) < 1e-3):
        return True
    else:
        return False
    
def blank_at_zero(tick,label):
    tolerance = 1.e-14
    if type(tick) is np.ndarray:
        index = np.where(abs(tick) < tolerance)
    else:
        temp = np.array(tick)
        index = np.where(abs(temp) < tolerance)
        del temp

    if np.size(index) == 0:
        raise ValueError('Array must span negative to positive values tick=',tick)
    else:
        index = index[0].item()
        label[index] = ''
        
def get_axis_tick_label(value):
    '''
    Get label for number on axis without trailing zeros.
    59.400000000000006
    will be returned as a string 
    '59.4'
    '''
    number_digits = 0
    if not use_sci_notation(value):
        label = str(value)
        
        # Get substring after period
        trailing = label.partition('.')[2]
        number_sigfig = 0
        if len(trailing) > 0:
            # Find number of non-zero digits after decimal
            number_sigfig = 1
            before = trailing[0]
            number_digits = 1
            go = True
            while go and number_digits < len(trailing):
                if trailing[number_digits] == before:
                    number_sigfig = number_digits - 1
                    if(number_sigfig > 5): go = False
                else:
                    before = trailing[number_digits]
                    number_sigfig = number_digits - 1
                number_digits+=1
    
        if number_digits == len(trailing): number_sigfig = number_digits

        # Round up the number to desired significant figures
        label = str(round(value, number_sigfig))
    else:
        label = "{:.1e}".format(value)

    return label

def pol2cart(phi, rho):
    '''
    make polar coords cartesian
    
    INPUTS:
    phi : polar angle counter-clockwise from x-axis in radians
    rho : radius
    
    OUTPUTS:
    x   : Cartesian x-coordinate
    y   : Cartesian y-coordinate
    '''

    x = np.multiply(rho, np.cos(phi))
    y = np.multiply(rho, np.sin(phi))
    return x, y