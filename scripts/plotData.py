import commands
import configIO
import os
import re
import numpy as np
from pylab import *
import itertools
import time

def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

def readColumnFromFile(col_idx, filename):
    f = open(filename)
    lines = f.readlines()

    if len(lines) == 0:
        return []

    # check if first line is a header
    first_idx = 0
    l = lines[0]
    l_s = l.split()
    if not isFloat(l_s[0]):
        #print l_s[0] + ' is not a digit'
        first_idx = 1

    objs = [0 for x in range(len(lines)-first_idx)]
    obj_idx = 0
    for li in range(first_idx, len(lines)):
        l = lines[li]
        l_s = l.split()
        if len(l_s) == 0:
            break
        if len(l_s) <= col_idx:
            obj = 0
        else:
            obj = float(l_s[col_idx])
        objs[obj_idx] = obj
        obj_idx = obj_idx + 1
    f.close()
    return objs

def getListCs():
    colors = ['b','g','r','c','m','y','k']
    styles = ['-','--','-.',':','.','o','v','^','v']
    listOfLists = [styles, colors]
    list_cs = list(itertools.product(*listOfLists))    
    return list_cs

def getLabels(dir_names, filename, max_delta_time = -1, dir_name_to_exclude = '', suffix = ''):

    current_time = time.time()

    filename_noext, ext = os.path.splitext(filename)
    filename_noext = filename_noext + suffix
    ldirs = os.listdir('.')
    ndirs = len(ldirs)

    labels = []

    for dir_name in dir_names:

        dir_name = re.escape(dir_name)

        for _dir in ldirs:

            if not os.path.isdir(_dir):
                continue
            m = re.search(dir_name, _dir)
            if m is None:
                continue

            if dir_name_to_exclude != '':
                m = re.search(dir_name_to_exclude, _dir)
                if m is not None:
                    continue

            fullpath = _dir + '/' + filename

            if not os.path.isfile(fullpath):
                continue

            if max_delta_time > 0:
                delta_time = current_time - os.path.getmtime(_dir)
                if delta_time > max_delta_time:
                    continue

            label_name = _dir
            labels.append(label_name)

    #common_label = os.path.commonprefix(labels)
    #for i in range(0,len(labels)):
    #    labels[i] = ''.join(labels[i].split(common_label))

    return labels


def plotData(dir_names, filename, max_delta_time = -1, dir_name_to_exclude = '', show_legend = 1, col_idx = 0, suffix = '', filename_multiCoeff = '', coeff_ylim = -1):

    current_time = time.time()

    filename_noext, ext = os.path.splitext(filename)
    filename_noext = filename_noext + suffix
    ldirs = os.listdir('.')
    ndirs = len(ldirs)

    fig_idx = 0
    list_objs = []
    labels = []

    for dir_name in dir_names:

        #dir_name = re.escape(dir_name)

        for _dir in ldirs:

            #_dir = re.escape(_dir)

            if not os.path.isdir(_dir):
                continue
            m = re.search(dir_name, _dir)
            if m is None:
                #print dir_name + ' not found in ' + _dir
                continue

            if dir_name_to_exclude != '':
                m = re.search(dir_name_to_exclude, _dir)
                if m is not None:
                    continue

            fullpath = _dir + '/' + filename

            if not os.path.isfile(fullpath):
                #print fullpath + ' is not a file'
                continue

            if os.path.getsize(fullpath) == 0:
                continue

            if max_delta_time > 0:
                delta_time = current_time - os.path.getmtime(_dir)
                if delta_time > max_delta_time:
                    continue

            label_name = _dir
            labels.append(label_name)

            if filename_multiCoeff != '':
                nConstraints = readColumnFromFile(2, _dir + '/' + filename_multiCoeff)
                objs0 = readColumnFromFile(col_idx, fullpath)
                #print 'nConstraints: ' + str(len(nConstraints)) + str(len(objs0)) + str(sum(nConstraints))
                nStepsPerOutput = len(nConstraints)/len(objs0)
                objs = [0 for x in range(int(sum(nConstraints)))] 
                idx_obj = 0
                idx_constraint = 0
                for i in range(0, len(objs0)):
                    for n in range(0,nStepsPerOutput):
                        for c in range(0, int(nConstraints[idx_constraint])):
                            objs[idx_obj] = objs0[i]
                            idx_obj = idx_obj + 1
                        idx_constraint = idx_constraint + 1

                # Check if using old file format
                if idx_obj == 0:
                    objs = readColumnFromFile(col_idx, fullpath)
            else:
                objs = readColumnFromFile(col_idx, fullpath)

            list_objs.append(objs)
            fig_idx = fig_idx + 1

    common_label = os.path.commonprefix(labels)

    for i in range(0,fig_idx):
        labels[i] = ''.join(labels[i].split(common_label))

    # Output curves
    if fig_idx > 0:
        list_cs = getListCs()
        figure(num=fig_idx, figsize=(12, 10), dpi=100)
        clf
        for i in range(0,len(list_objs)):
            l = len(list_objs[i])
            cs = list_cs[i][0] + list_cs[i][1]
            #y = list_objs[i][2:l]
            y = list_objs[i][0:l]
            plot(y, cs, linewidth=1.0)

            if coeff_ylim > 0:
                mean_y = np.mean(y)
                std_y = np.std(y)
                ylim([mean_y - coeff_ylim*std_y, mean_y + coeff_ylim*std_y])

        if show_legend == 1:
            legend(labels, loc='upper right')
        figtext(.3,.05,common_label)
        if col_idx == 0:
            output_file = os.path.basename(filename_noext) + '.png'
        else:
            output_file = os.path.basename(filename_noext) + str(col_idx) + '.png'
        print 'output_file ' + output_file
        try:
            savefig(output_file)
        except:
            pass
        clf
        close(fig_idx)

    return list_objs


def plotDataFromMultiColumns(dir_names, filename, max_delta_time = -1, dir_name_to_exclude = '', show_legend = 1, col_idx = [0], suffix = '', filename_multiCoeff = ''):

    current_time = time.time()

    filename_noext, ext = os.path.splitext(filename)
    filename_noext = filename_noext + suffix
    ldirs = os.listdir('.')
    ndirs = len(ldirs)

    fig_idx = 0
    list_objs = []
    labels = []

    for dir_name in dir_names:

        for _dir in ldirs:

            if not os.path.isdir(_dir):
                continue
            m = re.search(dir_name, _dir)
            if m is None:
                continue

            if dir_name_to_exclude != '':
                m = re.search(dir_name_to_exclude, _dir)
                if m is not None:
                    continue

            fullpath = _dir + '/' + filename

            if not os.path.isfile(fullpath):
                continue

            if max_delta_time > 0:
                delta_time = current_time - os.path.getmtime(_dir)
                if delta_time > max_delta_time:
                    continue

            label_name = _dir
            labels.append(label_name)

            col_idx = [1,4,5,6]

            if filename_multiCoeff != '':
                nConstraints = readColumnFromFile(2, _dir + '/' + filename_multiCoeff)
                objs_col1 = readColumnFromFile(col_idx[0], fullpath)
                objs_col4 = readColumnFromFile(col_idx[1], fullpath)
                objs_col5 = readColumnFromFile(col_idx[2], fullpath)
                objs_col6 = readColumnFromFile(col_idx[3], fullpath)
                objs0 = [0 for x in range(len(objs_col1))]
                for k in range(0, len(objs_col1)):
                    objs0[k] = (objs_col1[k] + objs_col4[k])/(objs_col5[k]+objs_col6[k])
                #print 'nConstraints: ' + str(len(nConstraints)) + str(len(objs0)) + str(sum(nConstraints))
                nStepsPerOutput = len(nConstraints)/len(objs0)
                objs = [0 for x in range(int(sum(nConstraints)))] 
                idx_obj = 0
                idx_constraint = 0
                for i in range(0, len(objs0)):
                    for n in range(0,nStepsPerOutput):
                        for c in range(0, int(nConstraints[idx_constraint])):
                            objs[idx_obj] = objs0[i]
                            idx_obj = idx_obj + 1
                        idx_constraint = idx_constraint + 1

                # Check if using old file format
                if idx_obj == 0:
                    objs_col1 = readColumnFromFile(col_idx[0], fullpath)
                    objs_col4 = readColumnFromFile(col_idx[1], fullpath)
                    objs_col5 = readColumnFromFile(col_idx[2], fullpath)
                    objs_col6 = readColumnFromFile(col_idx[3], fullpath)
                    objs = [0 for x in range(len(objs_col1))]
                    for k in range(0, len(objs_col1)):
                        objs[k] = (objs_col1[k] + objs_col4[k])/(objs_col5[k]+objs_col6[k])
            else:
                objs_col1 = readColumnFromFile(col_idx[0], fullpath)
                objs_col4 = readColumnFromFile(col_idx[1], fullpath)
                objs_col5 = readColumnFromFile(col_idx[2], fullpath)
                objs_col6 = readColumnFromFile(col_idx[3], fullpath)
                objs = [0 for x in range(len(objs_col1))]
                for k in range(0, len(objs_col1)):
                    objs[k] = (objs_col1[k] + objs_col4[k])/(objs_col5[k]+objs_col6[k])

            list_objs.append(objs)
            fig_idx = fig_idx + 1

    common_label = os.path.commonprefix(labels)

    for i in range(0,fig_idx):
        labels[i] = ''.join(labels[i].split(common_label))

    # Output curves
    if fig_idx > 0:
        list_cs = getListCs()
        figure(num=fig_idx, figsize=(12, 10), dpi=100)
        clf
        for i in range(0,len(list_objs)):
            l = len(list_objs[i])
            cs = list_cs[i][0] + list_cs[i][1]
            y = list_objs[i][2:l]
            plot(y, cs, linewidth=1.0)

        if show_legend == 1:
            legend(labels, loc='upper right')
        figtext(.3,.05,common_label)
        output_file = os.path.basename(filename_noext) + '.png'
        print 'output_file ' + output_file
        savefig(output_file)
        clf
        close(fig_idx)

    return list_objs

def plotWeights(dir_names, output_filename, max_delta_time = -1, dir_name_to_exclude = ''):

    current_time = time.time()

    ldirs = os.listdir('.')
    ndirs = len(ldirs)

    fig_idx = 0
    list_objs = []

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

            fullpath = _dir + '/parameter_vector0/'
            lfiles = os.listdir(fullpath)
            nfiles = len(lfiles)

            # Check if configuration file exists
            config_file_pattern = _dir[0:10]
            config_filename = ''
            for i in os.listdir(_dir):
                if i[0:len(config_file_pattern)] != config_file_pattern or i[len(i)-3:len(i)] != 'txt':
                    continue
                config_filename = _dir + '/' + i

            stepForParameterFiles = int(configIO.read(config_filename, 'stepForParameterFiles'))

            # plot last one only
            for i in range(nfiles-1, nfiles):
                fileId = i * stepForParameterFiles
                weight_filename = fullpath + 'iteration_' + str(fileId) + '.txt'
                #print 'weight ' + weight_filename
                if os.access(weight_filename, os.F_OK):
                    w = readColumnFromFile(0, weight_filename)
                    w = w[8:len(w)]

                    #fig = figure()
                    #fig = figure(num=None, figsize=(16, 12), dpi=80, facecolor='w', edgecolor='k')
                    fig = figure(num=None, figsize=(30, 14), dpi=160, facecolor='w', edgecolor=None)
                    ax = fig.add_subplot(1,1,1)
                    ind = range(0,len(w))
                    ax.bar(ind, w,align='center', edgecolor="none",linewidth=0.0)
                    #ax.bar(ind, w,align='center')
                    #ax.set_xticks(ind)
                    savefig(output_filename)
                    clf
