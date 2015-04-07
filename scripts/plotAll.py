#!/usr/bin/python

import matplotlib as mpl                                                                                                               
mpl.use('Agg')                                                                                                                 
import matplotlib.pyplot as plt

import os
import plotData
import plotCrossData
import plotLegend
import sys
from optparse import OptionParser
from pylab import *

def column(matrix, i):
    return [row[i] for row in matrix]

parser = OptionParser()
parser.add_option("-d", "--dir_name", dest="dir_name", default="training_files_", help="Directory name")
parser.add_option("-l", "--show_legend", dest="show_legend", default=1, help="Show legend")
parser.add_option("-n", "--dir_name_to_ignore", dest="dir_name_to_ignore", default="", help="Directory name to ignore")
parser.add_option("-o", "--output_dir", dest="output_dir", default="", help="Name of the directory where the html file will be copied")
parser.add_option("-t", "--max_delta_time", dest="max_delta_time", default=0.0, help="max_delta_time")
parser.add_option("-y", "--coeff_y_lim", dest="coeff_y_lim", default=-1.0, help="Coefficient for ylim")

(options, args) = parser.parse_args()

max_delta_time = 0
if options.max_delta_time != 0:
    max_delta_time = 60*24*int(options.max_delta_time)

coeff_y_lim = float(options.coeff_y_lim)

show_legend = int(options.show_legend)

output_dir = os.environ['HOME'] + '/public_html/' + options.output_dir

if not os.access(output_dir, os.F_OK):
    print 'Creating output directory ' + output_dir
    os.mkdir(output_dir)

dir_names = options.dir_name.split(',')

for d in xrange(0,len(dir_names)):
    dir_names[d] = dir_names[d].replace('+', '.*')

# name, col_idx, use_multi_coeff
var_names_with_attributes = [
    ['scores0/training_score', 9, 0],
    ['scores0/test_score', 9, 0],
    ['dscore', 0, 0],
    ['a_dscore', 0, 0],
    ['obj', 0, 0],
    ['m', 0, 0],
    ['learning_rate', 0, 0],
    ['norm_w', 0, 0],
    ['norm_dfy', 0, 0],
    ['loss', 0, 0],
    ['constraint_set_card', 0, 0],
    ['constraint_set_card', 1, 0],
    ['autostep_learning_rate', 0, 0],
    ['autostep_learning_rate_all_0',1, 0],
    ['autostep_learning_rate_all_1',1, 0],
    ['autostep_learning_rate_all_2',1, 0],
    ['autostep_linear_min',1, 0],
    ['autostep_quadratic_min',1, 0],
    ['d_slack',0, 0],
    ['max_slack_before',0, 0],
    ['fw_learning_rate',0, 0],
    ['fw_gap',0, 0],
    ['fw_dp',0, 0],
    ['fw_loss',0, 0],
    ['fw_norm_ws',0, 0],
    ['enforced_submodularity',0, 0],
    ['autostep_obj_0', 0, 0],
    ['autostep_obj', 0, 0],
    ['autostep_obj_plus', 0, 0],
    ['autostep_obj_minus', 0, 0],
    ['autostep_obj_qc', 0, 0],
    ['autostep_obj_qt', 0, 0],
    ['weights', 0, 0],
    ['legend', 0, 0]
    ];

var_names = column(var_names_with_attributes, 0)
col_idx = column(var_names_with_attributes, 1)
use_multi_coeff = column(var_names_with_attributes, 2)

weights_index = len(var_names) - 2
plotData.plotWeights(dir_names, 'weights.png', max_delta_time, options.dir_name_to_ignore)

legend_index = len(var_names) - 1
plotLegend.plotLegend(dir_names, 'legend.png', max_delta_time, options.dir_name_to_ignore)

filename_multiCoeff = 'constraint_set_card.txt'

labels =  plotData.getLabels(dir_names, var_names[0] + '.txt', max_delta_time, options.dir_name_to_ignore)

list_plots = []

# do not plot legend (need to use the plotLegend function called above)
for i in range(0,len(var_names) - 1):
    if use_multi_coeff[i]:
        l = plotData.plotData(dir_names, var_names[i] + '.txt', max_delta_time, options.dir_name_to_ignore, show_legend, col_idx[i],'',filename_multiCoeff, coeff_y_lim)
    else:
        l = plotData.plotData(dir_names, var_names[i] + '.txt', max_delta_time, options.dir_name_to_ignore, show_legend, col_idx[i],'', '', coeff_y_lim)
    list_plots.append(l)

# add empty plot for legend
list_plots.append([])

var_names_cross = [ 'd_norm_w',
                    'd_angle_w'
                    ];

#for i in range(0,len(var_names_cross)):
for i in range(0,0):
    l = plotCrossData.plotCrossData(dir_names, 'parameter_vector0/', i, var_names_cross[i] + '.png')
    list_plots.append(l)

###############
if 0:
    var_names_multi_cols = [
        'scores0/training_score',
        'scores0/test_score'
        ];

    for i in range(0,len(var_names_multi_cols)):
        l = plotData.plotDataFromMultiColumns(dir_names, var_names_multi_cols[i] + '.txt', max_delta_time, options.dir_name_to_ignore, show_legend, [1, 2],'',filename_multiCoeff)
###############

# Create figures with multiple variables
list_cs = plotData.getListCs()

# Pairs of plots to be plotted together on the same graph
pairs = [[0,1]]

pairs_add = []
#pairs_add = [[24,25], [24, 26], [24,25]]
operations = [1,1,0]

var_names_all = []
for i in range(0,len(var_names)):
    var_names_all.append(var_names[i])
for i in range(0,len(var_names_cross)):
    var_names_all.append(var_names_cross[i])

print 'var_names_all'
print var_names_all

# Display pairs of variables on the same graph
var_names_pairs = []
for p in range(0, len(pairs)):

    pair = pairs[p]
    i0 = pair[0]
    i1 = pair[1]

    figure(num=p, figsize=(12, 10), dpi=100)
    clf
    for ip in range(0,2):
        i = pair[ip]
        if (len(list_plots) > i) and len(list_plots[i]) > 1:
            l = len(list_plots[i][0])
            cs = list_cs[i][0] + list_cs[i][1]
            y = list_plots[i][0][2:l]
            plot(y, cs, linewidth=1.0)

    legend_name = os.path.basename(var_names_all[i0]) + '_' + os.path.basename(var_names_all[i1])
    var_names_pairs.append(legend_name)
    if show_legend == 1:
        legend(legend_name, loc='upper right')
    #figtext(.3,.05,common_label)
    output_file = legend_name + '.png'
    print 'output_file ' + output_file
    savefig(output_file)
    clf
    close(p)

for i in range(0,len(var_names_pairs)):
    var_names.append(var_names_pairs[i])
    var_names_all.append(var_names_pairs[i])

# Display addition of pairs of variables on the same graph
var_names_pairs_add = []
for p in range(0, len(pairs_add)):

    pair = pairs_add[p]

    figure(num=p, figsize=(12, 10), dpi=100)
    clf

    i0 = pair[0]
    i1 = pair[1]

    if (len(list_plots) <= i0) or (len(list_plots[i0]) == 0): 
        print 'Skipping pair ' + str(i0) + ',' + str(i1)
        continue

    cs = list_cs[i0][0] + list_cs[i0][1]
    y0 = list_plots[i0][0]
    y1 = list_plots[i1][0]

    ny = min(len(y0),len(y1))
    operation = operations[p]
    if operation == 0:
        y = [0 for i in xrange(ny)]
        for i in range(0,ny):
            y[i] = y0[i] - y1[i]
        y = [int(v>0) for v in y]
        #plot(y, cs, linewidth=1.0)
        plot(y, '.', linewidth=1.0)
    else:
        yt = [0 for i in xrange(ny)]
        for i in range(0,ny):
            yt[i] = y0[i] - y1[i]
        y = [0 for i in xrange(3)]
        y[0] = sum(v > 0 for v in yt)
        y[1] = sum(v == 0 for v in yt)
        y[2] = sum(v < 0 for v in yt)
        #print(y)
        fig = figure()
        ax = fig.add_subplot(1,1,1)
        ind = [1,2,3]
        ax.bar(ind, y,align='center')
        ax.set_xticks(ind)
        y_labels = ['>0','0','<0']
        ax.set_xticklabels(y_labels)

    legend_name = os.path.basename(var_names_all[i0]) + '_' + os.path.basename(var_names_all[i1]) + '_' + str(operation)
    var_names_pairs_add.append(legend_name)
    if show_legend == 1:
        legend(legend_name, loc='upper right')
    #figtext(.3,.05,common_label)
    output_file = legend_name + '.png'
    print 'output_file ' + output_file
    savefig(output_file)
    clf
    close(p)

for i in range(0,len(var_names_pairs_add)):
    var_names.append(var_names_pairs_add[i])
    var_names_all.append(var_names_pairs_add[i])

# Copy png file to output directory
if os.path.isdir(output_dir):
    cmd = 'cp '
    for i in range(0,len(var_names_all)):
        name = os.path.basename(var_names_all[i])
        if (i < len(col_idx)) and col_idx[i] != 0:
            name += str(col_idx[i])
        cmd = cmd + name + '.png '
    cmd = cmd + output_dir
    print cmd
    os.system(cmd)

f = open(output_dir + "plots.html", "w")
f.write('<HTML>\n')
f.write('<HEAD>\n')
f.write('<TITLE>Experiments - plots</TITLE>\n')
f.write('</HEAD>\n')
f.write('<BODY BGCOLOR="FFFFFF">\n')
f.write('<CENTER>\n')

f.write('<p>\n')

legend_html_text = ''
for j in range(0,len(var_names_all)):

    #if not (j == legend_index) and (j < len(list_plots)) and len(list_plots[j]) == 0:
    if not (j >= weights_index) and ((j < len(list_plots)) and len(list_plots[j]) == 0):
        #print 'Skip ' + var_names_all[j]
        continue

    name_j = os.path.basename(var_names_all[j])
    if (j < len(col_idx)) and col_idx[j] != 0:
            name_j += str(col_idx[j])

    legend_html_text = legend_html_text + '<a href="#' + name_j + '">' + name_j + '</a>&nbsp;&nbsp;&nbsp;'

f.write(legend_html_text)

for i in range(0,len(var_names)):

    #if not (i == legend_index) and (i < len(list_plots)) and len(list_plots[i]) == 0:
    if not (i >= weights_index) and ((i < len(list_plots)) and len(list_plots[i]) == 0):
        continue

    name = os.path.basename(var_names[i])
    if (i < len(col_idx)) and col_idx[i] != 0:
        name += str(col_idx[i])

    f.write('<p><a name="' + name + '">' + name + '</a></p>\n')
    f.write('<img src="' + name + '.png" align="bottom">\n')
    f.write('</br>\n')

    f.write(legend_html_text)

for i in range(0,len(var_names_cross)):

    name = os.path.basename(var_names_cross[i])

    f.write('<p><a name="' + name + '">' + name + '</a></p>\n')
    f.write('<img width=800 src="' + name + '.png" align="bottom">\n')
    f.write('</br>\n')

    f.write(legend_html_text)

f.write('</p>\n')
f.write('</CENTER>\n')
f.write('</BODY>\n')
f.write('</HTML>\n')
