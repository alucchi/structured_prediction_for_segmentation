import commands
import os
import re
from pylab import *
import itertools
import time

def plotLegend(dir_names, output_file, max_delta_time = -1, dir_name_to_exclude = ''):

    colors = ['b','g','r','c','m','y','k']
    styles = ['-','--','-.',':','.','o','v','^','v']
    listOfLists = [styles, colors]
    list_cs = list(itertools.product(*listOfLists))    

    current_time = time.time()

    ldirs = os.listdir('.')
    ndirs = len(ldirs)

    fig_idx = 0
    list_objs = []
    labels = []

    for dir_name in dir_names:

        for _dir in ldirs:

            if not os.path.isdir(_dir):
                continue
            #if _dir[0:len(dir_name)] != dir_name:
            #    continue
            m = re.search(dir_name, _dir)
            if m is None:
                continue
            #else:
            #    print 'Match: ' + _dir

            if dir_name_to_exclude != '':
                m = re.search(dir_name_to_exclude, _dir)
                if m is not None:
                    continue

            if max_delta_time > 0:
                delta_time = current_time - os.path.getmtime(_dir)
                if delta_time > max_delta_time:
                    continue

           #label_name = _dir[len(dir_name):len(_dir)-1]
            label_name = _dir
            labels.append(label_name)
            #labels[fig_idx] = label_name

            fig_idx = fig_idx + 1


    common_label = os.path.commonprefix(labels)

    #print 'common'
    #print labels
    print common_label

    for i in range(0,fig_idx):
        labels[i] = ''.join(labels[i].split(common_label))

    # Output curves
    if fig_idx > 0:
        figure(num=fig_idx, figsize=(12, 10), dpi=100)
        clf
        for i in range(0, len(labels)):
            cs = list_cs[i][0] + list_cs[i][1]
            plot(0, cs)

        legend(labels, loc='upper right')
        figtext(.3,.05,common_label)
        #axis([xmin, xmax, ymin, ymax])
        #axis([0, 8000, -1e12, 1e12])
        savefig(output_file)
        clf
        close(fig_idx)
