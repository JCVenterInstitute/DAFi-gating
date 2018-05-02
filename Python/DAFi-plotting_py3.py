#!/usr/bin/env python

import sys
import argparse
import glob
import errno
import math
import re
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.colors as cc
import matplotlib.lines as line2d
from matplotlib.ticker import MaxNLocator # added
import matplotlib.lines
import numpy as np
from color_palette import color_palette
from multiprocessing import Pool
import matplotlib.image as mpimg
from collections import defaultdict, OrderedDict
import traceback
import itertools
import warnings

_nsre = re.compile('([0-9]+)')
def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s)]   
#
# plotfig - function for plotting and writing each individual 2-D dot plot
#
def plotfig(fname, xdata, ydata, cdata, xlabel, ylabel, pops, lines, args):
    try:
        colorlist = args.colorlist
        sample_name = fname
                
        fig = plt.figure()
        plt.xlabel(xlabel, fontsize=args.sizeofaxislabels)
        plt.ylabel(ylabel, fontsize=args.sizeofaxislabels)
        ax = fig.add_subplot(111)
        c1 = cc.ColorConverter()
        c = c1.to_rgba("#D3D3D3",1.0)
        #c = c1.to_rgba("#000000",0.5)
        ax.patch.set_facecolor(c)
        ax.tick_params(axis='both',which='major',labelsize=10)
        ax.set_xlim([0,4096])
        ax.set_ylim([0,4096])
        ax.get_yaxis().set_tick_params(direction='in')
        ax.get_xaxis().set_tick_params(direction='in')
        ax.get_yaxis().set_ticks([0, 1024, 2048, 3072, 4096])
        ax.get_xaxis().set_ticks([0, 1024, 2048, 3072, 4096])
        ax.grid(args.gridshow, which='both')
        ax.set_axisbelow(True)
        if args.titleshow is True:
                        ax.set_title(sample_name, fontsize=args.sizeoftitle, color='red')
                        
        print "plotting figure..."
        plt.scatter(xdata,ydata,s=args.sizeofdot,alpha=args.alpha,c=cdata,edgecolors='none')
        
        if args.hidegatelines is False:
            for line in lines:
                ax.add_line(line2d.Line2D(line[0], line[1], linewidth=1, color=line[2]))
            
        if args.flocklegacy is False:
            poplist = ['Others'] + pops
            subcolorlist = colorlist[0:len(poplist)]
            pop_name=""
            
            if args.showparent:
                for i,pop in enumerate(pops):
                    if i>0: pop_name=pop_name+pop
            else:
                for i,pop in enumerate(pops):
                    pop_name=pop_name+pop
                    
            offset = 0
            for s,c in zip(poplist, subcolorlist):
                fig.text(0.92,0.3+offset," "+s+" ",color=c, fontsize=args.sizeofpoplabels)
                offset = offset + 0.02
        else:
            pop_name="FLOCK"
        
        F = plt.gcf()
        F.set_size_inches(6,6)
        F.set_dpi(300)
        header_names = xlabel+"_vs_"+ylabel
        png_file = sample_name + "_" + header_names + "_" + pop_name +".png"
        F.savefig(png_file,dpi=300)
        plt.close()
        return png_file
        
    except IOError as exc:
        if exc.errno != errno.EISDIR: # Do not fail if a directory is found, just ignore it.
            raise # Propagate other kinds of IOError.args
    except:
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))

#
#
#
def loadNp(openfile, nheaders, nlength):
    data = np.fromfile(openfile, sep=" ", dtype=int)
    #print data[0:nheaders]
    print np.prod(data.shape)
    data.shape = (nlength, nheaders)
    #print data[0,:]
    return data
        
#
# autoconfig_processfile - function to automatically parse configuration file and generate dot plot for each population gating
#
def autoconfig_processfile(name, pool_used, f_index, args):
    try:
        with open(name) as result_file: # No need to specify 'r': this is the default.as
            
            gates = args.config
            if(gates is not None):
                num_gates = len(gates)
            else:
                num_gates=0
                
            nameParts = name.split("/")
            
            if len(nameParts) > 1:
                originalName = nameParts[len(nameParts)-2]
            else:
                if args.name is not None:
                    originalName = args.name+str(f_index)
                else:
                    originalName = "NoNameSample"+str(f_index)
            print "Processing ", (originalName)
            
            colorlist = args.colorlist
            
            events = sum(1 for line in result_file) -1 #quickly determine number of events
            result_file.seek(0) #rewind to the beginning of file
            header = result_file.readline()
            header = header.strip()
            headers = header.split("\t") #parse the headers from the first line of input
            num_markers = len(headers) - 2
            
            # create a numpy array for faster data access
            if args.debug: print "Assigning data to numpy matrix"
            fcm = loadNp(result_file, len(headers), events)
            
            # find the start of pop info on fcs_results_all
            if args.flocklegacy is False:
                pop_offset=0
                for i,header in enumerate(headers):
                    if header == "pop1":
                        pop_offset=i-1
                if args.debug: print "Pop offset: ", pop_offset
            
            # parsing gates from configuration data
            if args.debug: print "Configuring gates from file"
            sub_results = []
            if pool_used == 0: inner_pool = Pool(processes=cores)
            for gate, config in gates.iteritems():
                xmarker=str(headers[config[1]-1])
                ymarker=str(headers[config[2]-1])
                startx=int((float(config[3])/200)*4096)
                starty=int((float(config[5])/200)*4096)
                endx=int((float(config[4])/200)*4096)
                endy=int((float(config[6])/200)*4096)
                parent_gate = int(config[7])
                
                lines=[]
                cluster_type=int(config[8])
                if cluster_type == 2:
                    print "slanted"
                else:
                    lines.append([(startx, endx), (starty, starty), colorlist[1]])
                    lines.append([(startx, startx), (starty, endy), colorlist[1]])
                    lines.append([(startx, endx), (endy, endy), colorlist[1]])
                    lines.append([(endx, endx), (starty, endy), colorlist[1]])
                    
                key=xmarker+"_vs_"+ymarker
                
                dim1 = xmarker
                dim2 = ymarker
                print dim1, dim2
                dim1_idx = 0
                dim2_idx = 0
                
                for i,marker in enumerate(headers):
                    if marker == dim1:
                        dim1_idx = i
                        print ("Feature 1: ", marker, i+1)
                    if marker == dim2:
                        dim2_idx = i
                        print ("Feature 2: ", marker, i+1)
                   
                header_names = key
                print header_names
                
                if args.flocklegacy is False:
                    poploc = int(config[0]) + pop_offset
                    parent_poploc = parent_gate + pop_offset
                    print "gate location: ", poploc
                    print "parent gate location: ", parent_poploc
                    if args.debug: print "iterate through events to find population members"
                    fcm[:,-1]=0
                    if args.showparent and (config[0] > 1):
                        fcm[:,-1]=(1-fcm[:,poploc])+(1-fcm[:,parent_poploc])
                    else:
                        fcm[:,-1]=(2-2*fcm[:,poploc])
                    
                if args.debug: print "sorting numpy array"
                if args.sort:
                    sfcm = fcm[np.argsort(fcm[:, -1])] #sort the data set based on population number
                else:
                    sfcm = fcm
                    
                if args.reversesort:
                    sfcm = sfcm[::-1]
                    
                #print sfcm[0, :]
                if args.debug: print "creating color array"
                cdata = []
                for a in sfcm[:, -1]:
                    cdata.append(colorlist[a])
                
                
                xdata = sfcm[:,dim1_idx]
                ydata = sfcm[:,dim2_idx]
                sample_name = originalName
                
                poplist = []
                if args.flocklegacy is False:
                    pop_name = gate
                    parent_pop_name = "pop"+str(parent_gate)
                    if args.showparent: poplist.append(parent_pop_name)
                    poplist.append(pop_name)
                else:
                    pop_name = "FLOCK"
                    
                
                png_file = sample_name + "_" + header_names + "_" + pop_name +".png"
                
                if pool_used == 0:
                    inner_pool.apply_async(plotfig, args=[sample_name, xdata, ydata, cdata, xmarker, ymarker, poplist, lines, args])
                else:
                    plotfig(sample_name, xdata, ydata, cdata, xmarker, ymarker, poplist, lines, args)
                if args.flocklegacy is False: 
                    sub_results.append([pop_name, [sample_name, png_file]])
                else:
                    sub_results.append([pop_name+str(gate), [sample_name, png_file]])
            if pool_used == 0:
                inner_pool.close()
                inner_pool.join()
            return sub_results
    except IOError as exc:
        if exc.errno != errno.EISDIR: # Do not fail if a directory is found, just ignore it.
            raise # Propagate other kinds of IOError.args
    except:
        #print >> sys.stderr, "Exception: %s" % str(e)
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))

def autoconfig_processfile2(name, pool_used, f_index, args):
    try:
        with open(name) as result_file: # No need to specify 'r': this is the default.as
            
            gates = args.config
            num_gates = len(gates)
            nameParts = name.split("/")
            
            if len(nameParts) > 1:
                originalName = nameParts[len(nameParts)-2]
            else:
                if args.name is not None:
                    originalName = args.name+str(f_index)
                else:
                    originalName = "NoNameSample"+str(f_index)
            print "Processing ", (originalName)
            
            colorlist = args.colorlist
            
            events = sum(1 for line in result_file) -1 #quickly determine number of events
            result_file.seek(0) #rewind to the beginning of file
            header = result_file.readline()
            header = header.strip()
            headers = header.split("\t") #parse the headers from the first line of input
            headers = filter(None, headers)
	    num_markers = len(headers) - 2
            
            # create a numpy array for faster data access
            if args.debug: 
		print "Assigning data to numpy matrix"
		print "header length: ", len(headers)
		print headers
		print "events count: ", events
            fcm = loadNp(result_file, len(headers), events)
            
            # find the start of pop info on fcs_results_all
            if args.flocklegacy is False:
                pop_offset=0
                for i,header in enumerate(headers):
                    if header == "pop1":
                        pop_offset=i-1
                if args.debug: print "Pop offset: ", pop_offset
            
            axis_popIndexDict = defaultdict(list)
            
            print "Configuring axises from gate configuration file"
            axises=[]
            composite_axis=0
            last_xmarker=""
            last_ymarker=""
            last_parent=0
            for pop, config in gates.items():
                #pop="pop"+str(i+1)
                #config=gates.get(pop)
                
                xmarker=str(headers[config[1]-1])
                ymarker=str(headers[config[2]-1])
                startx=int((float(config[3])/200)*4096)
                starty=int((float(config[5])/200)*4096)
                endx=int((float(config[4])/200)*4096)
                endy=int((float(config[6])/200)*4096)
                parent=int(config[7])
                key="axis"+str(composite_axis)
                if (xmarker != last_xmarker) or (ymarker != last_ymarker) or (parent != last_parent):
                    composite_axis=composite_axis+1
                    key="axis"+str(composite_axis)
                    axises.append([xmarker,ymarker,key])
                axis_popIndexDict[key].append(pop)
                last_xmarker=xmarker
                last_ymarker=ymarker
                last_parent=parent
                
            num_axises = len(axises)
            #print "axis_popIndexDict: ", axis_popIndexDict
            #print axises
            
            cols = int(min(num_axises,4))
            rows = int(max(int(math.ceil(num_axises / float(cols))),1))
            
            print "Iterating through feature pairs"
            sub_results = []
            composite_list = []
            if pool_used == 0: inner_pool = Pool(processes=cores)
            for w,mpair in enumerate(axises):
                lines = []
                dim1 = mpair[0]
                dim2 = mpair[1]
                #print dim1, dim2
                dim1_idx = 0
                dim2_idx = 0
                for i,marker in enumerate(headers):
                    if marker == dim1:
                        dim1_idx = i
                        print ("Feature 1: ", marker, i+1)
                    if marker == dim2:
                        dim2_idx = i
                        print ("Feature 2: ", marker, i+1)
                        
                header_names = dim1 + "_vs_" + dim2
                axis_name = mpair[2]
                print header_names, axis_name
                
                if args.flocklegacy is False:
                    fcm[:,-1]=0 #reset color mapping in np matrix
                    pops = axis_popIndexDict.get(axis_name)
                    if args.debug: print "iterate through events to find population members"
                    poplist_colors = []
                    for i, pop in enumerate(pops):
                        config=gates.get(pop)
                        found_pop = config[0]
                        xmarker=str(headers[config[1]-1])
                        ymarker=str(headers[config[2]-1])
                        startx=int((float(config[3])/200)*4096)
                        starty=int((float(config[5])/200)*4096)
                        endx=int((float(config[4])/200)*4096)
                        endy=int((float(config[6])/200)*4096)
                        if i==0: parent_gate = int(config[7])
                        cluster_type=int(config[8])
                        pop_loc = found_pop + pop_offset
                        parent_poploc = parent_gate + pop_offset
                        
                        
                        if (xmarker==dim1) and (ymarker==dim2):
                            if cluster_type == 2:
                                print "slanted"
                            else:
                                lines.append([(startx, endx), (starty, starty), colorlist[i+2]])
                                lines.append([(startx, startx), (starty, endy), colorlist[i+2]])
                                lines.append([(startx, endx), (endy, endy), colorlist[i+2]])
                                lines.append([(endx, endx), (starty, endy), colorlist[i+2]])
                        elif (xmarker==dim2) and (ymarker==dim1):
                            if cluster_type == 2:
                                print "slanted"
                            else:
                                lines.append([(starty, endy), (startx, startx), colorlist[i+2]])
                                lines.append([(starty, starty), (startx, endx), colorlist[i+2]])
                                lines.append([(starty, endy), (endx, endx), colorlist[i+2]])
                                lines.append([(endy, endy), (startx, endx), colorlist[i+2]])
                        
                        if args.showparent and (config[0] > 1) and i == 0:
                            fcm[:,-1]=(1-fcm[:,pop_loc])+(1-fcm[:,parent_poploc])
                        else:
                            fcm[:,-1]=np.maximum(fcm[:,-1], (1-fcm[:,pop_loc])*(i+2))
                                
                if args.debug: print "sorting numpy array"
                if args.sort:
                    sfcm = fcm[np.argsort(fcm[:, -1])] #sort the data set based on population number
                else:
                    sfcm = fcm
                    
                if args.reversesort:
                    sfcm = sfcm[::-1]
                    
                print "creating color array"
                cdata = []
		for a in sfcm[:, -1]:
			try:	
	                	cdata.append(colorlist[a])
			except Exception, err:
				print "ERROR in index: ", a
				sys.stderr.write('Error: %sn' % str(err))
				return 1   

                xdata = sfcm[:,dim1_idx]
                ydata = sfcm[:,dim2_idx]
                sample_name = originalName
                
                poplist = []
                if args.flocklegacy is False:
                    parent_pop_name = "pop"+str(parent_gate)
                    if args.showparent: 
                        poplist.append(parent_pop_name)
                    else:
                        poplist.append("")
                    poplist = poplist + pops
                    
                    pop_name=""
                    if args.showparent:
                        for i,pop in enumerate(poplist):
                            if i>0: pop_name=pop_name+pop
                    else:
                        for i,pop in enumerate(poplist):
                            pop_name=pop_name+pop
                else:
                    pop_name="FLOCK"
                        
                png_file = sample_name + "_" + header_names + "_" + pop_name +".png"
                if pool_used == 0:
                    inner_pool.apply_async(plotfig, args=[sample_name, xdata, ydata, cdata, dim1, dim2, poplist, lines, args])
                else:
                    png_file = plotfig(sample_name, xdata, ydata, cdata, dim1, dim2, poplist, lines, args)
                sub_results.append([axis_name, [sample_name, png_file]])
                composite_list.append([axis_name, png_file])
            if pool_used == 0:
                inner_pool.close()
                inner_pool.join()
            if args.gatescomposite:
                print "Number of rows and columns for the composite plots (gates): ", rows, cols    
                compose_from_fig(composite_list, rows, cols, num_axises, sample_name, " ", " ", 0, args)
            
            return sub_results
    except IOError as exc:
        if exc.errno != errno.EISDIR: # Do not fail if a directory is found, just ignore it.
            raise # Propagate other kinds of IOError.args
    except Exception,e:
        print >> sys.stderr, "Exception: %s" % str(e)
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))

#
# processfile - default function to process command line specifications of populations and markers to be plotted
#
def processfile(name, pool_used, f_index, args):
    try:
        with open(name) as result_file: # No need to specify 'r': this is the default.as
            gates = args.config
            if(gates is not None): num_gates = len(gates)
            nameParts = name.split("/")
            if args.name is not None:
                originalName = args.name+str(f_index)
            elif (args.flocklegacy is True):
                originalName = nameParts[len(nameParts)-1]
            else:
                if len(nameParts) > 1:
                    originalName = nameParts[len(nameParts)-2]
                else:
                    originalName = "NoNameSample"+str(f_index)
            print "Processing ", (originalName)
            
            # retrieve command line options from args
            colorlist = args.colorlist
            mpairs = args.markers
            
            events = sum(1 for line in result_file) -1 #quickly determine number of events
            
            result_file.seek(0) #rewind file to beginning
            header = result_file.readline()
            header = header.strip()
            headers = header.split("\t")
            markers = headers
            num_markers = len(headers) - 2
            
            if ((args.markers is None) and (args.flocklegacy is True)):
                mpairs = list(itertools.permutations(markers[:-2], 2))
            
            
            num_pairs = len(mpairs)
            
            #
            # create a numpy array
            #
            fcm = loadNp(result_file, len(headers), events)
            
            #
            # Check legacy FLOCK or new DAFI results
            #
            pops = args.listofpop
            if args.flocklegacy is False:
                pop2DList = [[i, j] for i, j in enumerate(pops)]
                # locate position of the specified populations in the input file columns
                for i,marker in enumerate(markers):
                    for j,k in enumerate(pops):
                        if marker == k:
                            pop2DList[j][0] = i
                # create a dictionary for quick population order lookup
                plist = [int(x[0]) for x in pop2DList]
                pop2DIndex = set(plist)
                pdict = {e:(plist.index(e)+1) for e in plist}
                pdict.update({0:0})
                
                # find the start of pop info on fcs_results_all
                pop_offset=0
                for i,header in enumerate(headers):
                    if header == "pop1":
                        pop_offset=i-1
                if args.debug: print "Pop offset: ", pop_offset
                
                if args.debug: print ("Events: ",events)
                if args.debug: print ("Number of Markers: ",num_markers)
                
                fcm[:,-1]=0
                for i, popIndex in enumerate(pop2DIndex):
                    print i, popIndex
                    fcm[:,-1] = np.maximum(fcm[:,-1], (1-fcm[:,popIndex])*(popIndex))
            
            # Sort the matrix data based on the population column (last column)
            if args.sort:
                sfcm = fcm[np.argsort(fcm[:, -1])] #sort the data set based on population number
            else:
                sfcm = fcm
                
            if args.reversesort:
                sfcm = sfcm[::-1]
            
            #generate color array from sorted numpy matrix
            colors = []
            for a in sfcm[:, -1]:
                if args.flocklegacy:
                    colors.append(colorlist[a])
                else:
                    colors.append(colorlist[pdict[a]])
            
            #calculate cols and rows for composite plot
            cols = int(min(num_pairs,4))
            rows = int(max(int(math.ceil(num_pairs / float(cols))),1))
            
            sub_results = []
            composite_list = []
            #iterate through all 2D dot plot pairs specified
            for w,mpair in enumerate(mpairs):
                
                dim1 = mpair[0]
                dim2 = mpair[1]
                if args.debug: print dim1, dim2
                dim1_idx = 0
                dim2_idx = 0
                
                lines=[]
                if((gates is not None) and (pops is not None)):
                    for i, pop in enumerate(pops):
                        config=gates.get(pop)
                        xmarker=str(headers[config[1]-1])
                        ymarker=str(headers[config[2]-1])
                        startx=int((float(config[3])/200)*4096)
                        starty=int((float(config[5])/200)*4096)
                        endx=int((float(config[4])/200)*4096)
                        endy=int((float(config[6])/200)*4096)
                        parent_gate = int(config[7])
                        cluster_type=int(config[8])
                        if (xmarker==dim1) and (ymarker==dim2):
                            if cluster_type == 2:
                                print "slanted"
                            else:
                                lines.append([(startx, endx), (starty, starty), colorlist[pdict[config[0]+pop_offset]]])
                                lines.append([(startx, startx), (starty, endy), colorlist[pdict[config[0]+pop_offset]]])
                                lines.append([(startx, endx), (endy, endy), colorlist[pdict[config[0]+pop_offset]]])
                                lines.append([(endx, endx), (starty, endy), colorlist[pdict[config[0]+pop_offset]]])
                        elif (xmarker==dim2) and (ymarker==dim1):
                            if cluster_type == 2:
                                print "slanted"
                            else:
                                lines.append([(starty, endy), (startx, startx), colorlist[pdict[config[0]+pop_offset]]])
                                lines.append([(starty, starty), (startx, endx), colorlist[pdict[config[0]+pop_offset]]])
                                lines.append([(starty, endy), (endx, endx), colorlist[pdict[config[0]+pop_offset]]])
                                lines.append([(endy, endy), (startx, endx), colorlist[pdict[config[0]+pop_offset]]])
                                
                for i,marker in enumerate(markers):
                    if marker == dim1:
                        dim1_idx = i
                        print ("Feature 1: ", marker, i)
                    if marker == dim2:
                        dim2_idx = i
                        print ("Feature 2: ", marker, i)
                        
                xdata = sfcm[:,dim1_idx]
                ydata = sfcm[:,dim2_idx]
                cdata = colors
                header_names = dim1 + "_vs_" + dim2
                sample_name = originalName
                
                if args.flocklegacy is False:
                    pop_name=""
                    if args.showparent:
                        for i,pop in enumerate(pops):
                            if i>0: pop_name=pop_name+pop
                    else:
                        for i,pop in enumerate(pops):
                            pop_name=pop_name+pop
                
                png_file = plotfig(sample_name, xdata, ydata, cdata, dim1, dim2, pops, lines, args)
                print png_file
                sub_results.append([header_names, [sample_name, png_file]])
                composite_list.append([header_names, png_file])
            
            if args.pairscomposite:
                print "Number of rows and columns for the composite plots (markers): ", rows, cols
                print composite_list
                compose_from_fig(composite_list, rows, cols, num_pairs, sample_name, " ", " ", 0, args)
            return sub_results
    except IOError as exc:
        if exc.errno != errno.EISDIR: # Do not fail if a directory is found, just ignore it.
            raise # Propagate other kinds of IOError.args
    except:
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))

#
# compose_from_fig - function to generate composite plot from figures pre-generated via individual plotfig calls
#
def compose_from_fig(data, rows, cols, num_images, title, xlabel, ylabel, spaceing, args):
    try:
        
        #for autoparsing of configuration file, population legends are plotted at individual figure level
        if args.listofpop is not None: 
            pops = args.listofpop
        else:
            pops=[]
        
        if args.flocklegacy is False:        
            pop_name=""
            if args.showparent:
                for i,pop in enumerate(pops):
                    if i>0: pop_name=pop_name+pop
            else:
                for i,pop in enumerate(pops):
                    pop_name=pop_name+pop
        else:
            pop_name="FLOCK"
                
        colorlist = args.colorlist
        dpi = args.dpiofoutput
        currentRows = rows
        processingRows = 0
        processedRows = 0
        pageCount = 0
        startOffset = 0
        
        print "Number of images to compose: "+str(num_images)+", number of rows: "+str(rows)
        while (currentRows > 0):
            
            if currentRows>4:
                processingRows=4
            else:
                processingRows=currentRows
            
            if args.debug: print "currentRows: "+str(currentRows)+", processedRows: "+str(processedRows)+", processing "+str(processingRows)+" of rows, on page "+str(pageCount)
                
            f, axarr = plt.subplots(processingRows, cols)
            f.set_size_inches(8, 2+(2*processingRows))
            #plt.tight_layout(w_pad=2.0, h_pad=1.0)
            f.set_dpi(dpi)
            f.suptitle(title, size=20)
            
            
            c1 = cc.ColorConverter()
            c = c1.to_rgba("#D3D3D3",1.0)
            
            for i, ax in enumerate(axarr.flat, start=0):
                rowIndex = i/cols
                colIndex = i%cols
                imageIndex=(i+startOffset)
                if args.debug: print "Processing image: ", rowIndex, colIndex, imageIndex
                ax.set_axis_off()
                
                if imageIndex < num_images:
                    returnObject = data[imageIndex]
                    sname = returnObject[0]
                    fname = returnObject[1]
                    #snames = sname.split("_vs_")
                    img=mpimg.imread(fname)
                    ax.imshow(img)
                    if args.titleshow is True:
                        ax.set_title(sname, fontsize=6, color='red')
                
            f.subplots_adjust(wspace=spaceing)
            f.subplots_adjust(hspace=spaceing)
                
            png_file = title+pop_name+"_composite_"+str(pageCount)+".png"
            f.text(0.5, 0.04, xlabel, ha='center') #x-axis label
            f.text(0.02, 0.5, ylabel, va='center', rotation='vertical') #y-axis label
            
            if args.flocklegacy is False:
                if len(pops) > 0:
                    poplist = ['Others'] + pops
                    subcolorlist = colorlist[0:len(poplist)]
                    
                    offset = 0
                    for s,c in zip(poplist, subcolorlist):
                        f.text(0.92,0.3+offset," "+s+" ",color=c, fontsize=5)
                        offset = offset + 0.02
                    
            if args.debug: print (png_file)
            f.savefig(png_file,dpi=1200)
            
            processedRows=processedRows+processingRows
            currentRows=currentRows-processingRows
            startOffset = processedRows*cols
            pageCount=pageCount+1
            plt.close()
    except Exception, e:
        print >> sys.stderr, "Exception: %s" % str(e)
        print "Problem file: "+title
        print "Length of data package: ", len(data)
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))

#
# sample_composite - function to compose the pre-generated figures for each sample on a grid
#
def sample_composite(convResults, num_sources, p_args):
    dataDict = defaultdict(list)
            
    for i in range(len(convResults)):
        for j in range(len(convResults[i])):
            data_name = convResults[i][j][0]
            if p_args.debug: print "Processing data pack ", data_name
            dataDict[data_name].append(convResults[i][j][1])
            #templist = dataDict.get(data_name)
            # for k,data_test in enumerate(templist):
            #     print "added ", data_test[0], " to ", data_name
    scomp_pool = Pool(processes=cores)
    try:
        if num_sources > 1:
            
            for w,currentPair in enumerate(dataDict):
                samplelist = dataDict.get(currentPair)
                sorted_results = sorted(samplelist,key=lambda x: (x[0]))
                
                num_samples = len(samplelist)
                cols = int(min(num_samples,4))
                rows = int(max(int(math.ceil(num_samples / float(cols))),1))
                print "Number of rows and columns for the composite plots (samples): ", rows, cols
    
                dim1=""
                dim2=""
                headerstr = str(currentPair).split("_vs_")
                if len(headerstr) > 1:
                    dim1 = str(headerstr[0])
                    dim2 = str(headerstr[1])
                
                print "Processing composite for "+currentPair
                scomp_pool.apply_async(compose_from_fig, args=[sorted_results, rows, cols, num_samples, currentPair, dim1, dim2, 0, p_args])
            scomp_pool.close()
            scomp_pool.join()
        else:
            #procedure for composing all populations of a single sample
            composite_list = []
            sample_name = ""
            dictKeys = dataDict.keys()
            dictKeys.sort(key=natural_sort_key)
            print "sorted keys: ", dictKeys
            for w,key in enumerate(dictKeys):
                samplelist = dataDict.get(key)
                if w==0: sample_name=samplelist[0][0]
                samplelist[0][0]=key #replace file name with population name
                composite_list=composite_list + samplelist
            #sorted_results = sorted(composite_list,key=lambda x: (x[0]))
            num_figs = len(composite_list)
            cols = int(min(num_figs,4))
            rows = int(max(int(math.ceil(num_figs / float(cols))),1))
            print "Number of rows and columns for the composite plots (samples): ", rows, cols
            compose_from_fig(composite_list, rows, cols, num_figs, sample_name, "", "", 0, p_args)
        plt.close("all")
    except KeyboardInterrupt:
        scomp_pool.terminate()
        scomp_pool.join()
        sys.exit(1)
    except Exception, e:
        scomp_pool.terminate()
        scomp_pool.join()
        print >> sys.stderr, "Exception: %s" % str(e)
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))
    
#
# marker_pair - function for parsion customized marker pair type for argparse
#
def marker_pair(s):
    try:
        xmarker, ymarker = map(str, s.split(','))
        return [xmarker, ymarker]
    except:
        raise argparse.ArgumentTypeError("marker pair must be xmarker,ymarker")
#
# config_objects - function for parsion customized configuration file type for argparse
#        
def config_objects(s):
    try:
        with open(s) as config_file:
            config_file.seek(0)
            gates=OrderedDict()
            for line in config_file:
                line = line.strip()
                gate = line.split("\t")
                gates.update({"pop"+str(gate[0]):[int(gate[0]), int(gate[1]), int(gate[2]), int(gate[3]), int(gate[4]), int(gate[5]), int(gate[6]), int(gate[7]), int(gate[8])]})
            return gates
    except:
        raise argparse.ArgumentTypeError("Error parsing configuration file")

#
# main function - argparse and directing multithreading
#
if __name__ == '__main__':
    #argparse for parsing command line
    parser = argparse.ArgumentParser(description='A flow cytometry 2D dot-plot generator from DAFi/FLOCK computational cell filtering results.')
    parser.add_argument('-f','--fileinput', help='<Required>Input DAFi/FLOCK result file name (use wildcard to do batch processing', required=True)
    parser.add_argument('--config',help='[Optional]Config file for automated marker/pop selection', default=None, required=False, type=config_objects)
    parser.add_argument('-m','--markers',help='<Required>Specify pair(s) of markers to be plotted. (e.g. CD3,CD4 CD4,Tet)', default=None, action='append', type=marker_pair, required=False)
    parser.add_argument('-l','--listofpop',help='<Required>List of cell populations to be highlighted together (e.g. "-l pop5 pop10")', default=None, nargs='+', required=False, type=str)
    parser.add_argument('--poplist',help='<Required>List of cell populations to be plotted separately (e.g. "--poplist pop5 pop10")', default=None, nargs='+', required=False, type=str)
    parser.add_argument('-c','--colorlist',help='Set color list as palette (please provide a list as long as the list of population to be displayed', default=color_palette, nargs="+", required=False, type=str)
    parser.add_argument('-p','--pthread',help='Number of available cores (default: 2)', default=2, required=False, type=int)
    parser.add_argument('-t','--titleshow',help='Display file title for each subplot', action="store_true", default=False, required=False)
    parser.add_argument('--legendshow',help='Display color legend', action="store_true", default=False, required=False)
    parser.add_argument('-s','--sizeofdot',help='Set point size for each plotted dot', default=1.0, required=False, type=float)
    parser.add_argument('--sizeofaxislabels',help='Set size for axis text labels', default=10, required=False, type=int)
    parser.add_argument('--sizeoftitle',help='Set font size for title', default=10, required=False, type=int)
    parser.add_argument('--sizeofpoplabels',help='Set size for population labels', default=8, required=False, type=int)
    parser.add_argument('-d','--dpiofoutput',help='Set number of dots per inch resolution for output images', default=300, required=False, type=int)
    parser.add_argument('-a','--alpha',help='Set alpha blending value of marker', default=1.0, required=False, type=float)
    parser.add_argument('--gridshow',help='Display grid for each plot', action="store_true", default=False, required=False)
    parser.add_argument('--samplescomposite',help='Generate composite plot across all files', action="store_true", default=False, required=False)
    parser.add_argument('--pairscomposite',help='Generate composite plot across all marker pairs', action="store_true", default=False, required=False)
    parser.add_argument('--gatescomposite',help='Generate composite plot across all marker pairs', action="store_true", default=False, required=False)
    parser.add_argument('--individual',help='Generate individual plots of marker pair', action="store_true", default=False, required=False)
    parser.add_argument('--debug',help='debug info output', action="store_true", default=False, required=False)
    parser.add_argument('--showparent',help='Show parent population', action="store_true", default=False, required=False)
    parser.add_argument('--showmultigates',help='Show multiple gates on the plot', action="store_true", default=False, required=False)
    parser.add_argument('--hidegatelines',help='Hide lines of the gates boundaries', action="store_true", default=False, required=False)
    parser.add_argument('--sort',help='enable sorting by populations', action="store_true", default=False, required=False)
    parser.add_argument('--reversesort',help='reverse sorting by populations', action="store_true", default=False, required=False)
    parser.add_argument('--flocklegacy',help='support for original FLOCK output', action="store_true", default=False, required=False)
    parser.add_argument('--name',help='Manual naming of sample', default=None, required=False, type=str)
    
    p_args = parser.parse_args()
    
    path = p_args.fileinput
    warnings.filterwarnings("ignore")    
    # check if either configuration file or manual configuration is specified
    if (p_args.flocklegacy is False):
        if (((p_args.config is None) and (p_args.markers is None)) or ((p_args.config is None) and (p_args.listofpop is None) and (p_args.markers is None))):
            print >> sys.stderr, "Must specify either the configuration file or manual list of populations and marker pairs (see help)"
            sys.exit(1)
    else:
        if ((p_args.config is not None) or (p_args.listofpop is not None)):
            print >> sys.stderr, "Can't specify configuration and/or populations while in legacy FLOCK mode (see help)"
            sys.exit(1)
    
    pops = p_args.listofpop
    colorlist = p_args.colorlist
    cores = p_args.pthread
    sizeP = p_args.sizeofdot
    dpi = p_args.dpiofoutput
    
    #parse out wildcard selected files
    files = glob.glob(path)
    num_files = len(files)
    
    if num_files > 1:
        file_pool = Pool(processes=cores)
        try:
            results = []
            for i, name in enumerate(files): # 'file' is a builtin type, 'name' is a less-ambiguous variable name.
                if ((p_args.markers is None) and (p_args.flocklegacy is False)):
                    if p_args.showmultigates:
                        proc = file_pool.apply_async(autoconfig_processfile2, args=[name, min(24, num_files), i, p_args])
                    else:
                        proc = file_pool.apply_async(autoconfig_processfile, args=[name, min(24, num_files), i, p_args])
                else:
                    proc = file_pool.apply_async(processfile, args=[name, min(24, num_files), i, p_args])
                results.append(proc)
            file_pool.close()
            file_pool.join()
        except KeyboardInterrupt:
            file_pool.terminate()
            file_pool.join()
            sys.exit(1)
        except Exception, e:
            file_pool.terminate()
            file_pool.join()
            print >> sys.stderr, "Exception: %s" % str(e)
            raise Exception("".join(traceback.format_exception(*sys.exc_info())))
            
        plt.close("all")
        
        #generating composite at samples level
        if p_args.samplescomposite:
            
            convResults = []
            for result in results:
                tempResult = result.get()
                if p_args.debug: print tempResult[0][0]
                convResults.append(tempResult)
            
            sample_composite(convResults, num_files, p_args)
	sys.exit(0)
    else:
        try:
            print "Processing only 1 file........................."
            if (((p_args.listofpop is None) and (p_args.flocklegacy is False)) or ((p_args.markers is None) and (p_args.flocklegacy is False))):
                if p_args.showmultigates:
                    convResults = [autoconfig_processfile2(files[0], 0, 0, p_args)]
                else:
                    convResults = [autoconfig_processfile(files[0], 0, 0, p_args)]
            else:
                convResults = [processfile(files[0], 0, 0, p_args)]
                
            if p_args.samplescomposite:
                sample_composite(convResults, 0, p_args)
            sys.exit(0)    
        except KeyboardInterrupt:
            sys.exit(1)
        except Exception, e:
            print >> sys.stderr, "Exception: %s" % str(e)
            raise Exception("".join(traceback.format_exception(*sys.exc_info())))
	sys.exit(0)
sys.exit(0)
