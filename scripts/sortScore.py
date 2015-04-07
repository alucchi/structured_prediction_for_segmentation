#!/usr/bin/python

import commands
import configIO
import itertools
import math
import os
import re
from pylab import *
import time
from optparse import OptionParser

def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

def trunc(f, n):
    '''Truncates/pads a float f to n decimal places without rounding'''
    slen = len('%.*f' % (n, f))
    return float(str(f)[:slen])

#print 'Usage: plotScore [dir_name] [training=0,test=1] [EM_3d=0,EM_2d=1,MSRC=2]'

parser = OptionParser()
parser.add_option("-d", "--dir_name", dest="dir_name", default="training_files_", help="Directory name")
parser.add_option("-i", "--name_to_ignore", dest="name_to_ignore", default="", help="Name to ignore")
parser.add_option("-m", "--max_iteration_to_show", dest="max_iteration_to_show", default=-1, help="Maximum iterations to show")
parser.add_option("-r", "--round_max_score", dest="round_max_score", default=0, help="round_max_score")
parser.add_option("-s", "--score_type", dest="score_type", default=0, help="Score type. 0 = training, 1 = test, 2 = test with best training score")
parser.add_option("-t", "--max_delta_time", dest="max_delta_time", default=-1, help="max_delta_time")
parser.add_option("-x", "--sorting_index", dest="sorting_index", default=-1, help="Index to use for sorting")

(options, args) = parser.parse_args()

dir_name = options.dir_name
dir_name = dir_name.replace('+', '.*')

max_iteration_to_show = int(options.max_iteration_to_show)
sorting_index = int(options.sorting_index)
max_delta_time = int(options.max_delta_time)
showRecentFilesOnly = max_delta_time != -1
name_to_ignore = options.name_to_ignore
round_max_score = int(options.round_max_score) == 1

score_filename = 'training_score.txt'
output_png_file = 'training_scores.png'
score_type = int(options.score_type)
if score_type == 1:
    #print 'display test scores'
    score_filename = 'test_score.txt'
    output_png_file = 'test_scores.png'
#else:
#    print 'display training scores'

showRunningTime = False

fig_idx = 1
list_max_scores = []
list_iteration_ids = []
list_n_iterations = []
list_times = []
labels = []
list_step_ids = []

list_files = []
recursive_mode = False
if recursive_mode == True:
    for _dir in os.listdir('.'):
        if os.path.isdir(_dir):
            for sub_dir in os.listdir(_dir):
                list_files.append(_dir + '/' + sub_dir)
        else:
            list_files.append(_dir)
else:
    for _dir in os.listdir('.'):
        list_files.append(_dir)
        
for _dir in list_files:
    if not os.path.isdir(_dir):
        continue
    m = re.search(dir_name, _dir)
    if m is None:
        continue
    #else:
    #    print 'Match: ' + _dir
    if name_to_ignore != '':
        m = re.search(name_to_ignore, _dir)
        if m is not None:
            continue

    if showRecentFilesOnly:
        delta_time = time.time() - os.path.getmtime(_dir)
        if delta_time > max_delta_time:
            continue

    output_dir = _dir + '/scores0/'
    if not os.path.isdir(output_dir):
        continue

    max_score = 0
    filename_max_score = ''
    scores = [0]
    idx_score = 0

    # Check if configuration file exists
    config_file_pattern = _dir[0:10]
    config_filename = ''
    for i in os.listdir(_dir):
        if i[0:len(config_file_pattern)] != config_file_pattern or i[len(i)-3:len(i)] != 'txt':
            continue
        config_filename = _dir + '/' + i

    stepForOutputFiles = 1
    if config_filename != '':
        stepForOutputFiles_str = configIO.read(config_filename, 'stepForOutputFiles')
        stepForOutputFiles = int(stepForOutputFiles_str)

    # Check if score file exists 
    filename = ''
    for i in os.listdir(output_dir):
        if i[0:len(score_filename)] != score_filename:
            continue
        filename = output_dir + i

    if filename != '':
        f = open(filename)
        lines = f.readlines()

        if sorting_index == -1:
            if _dir.find('MSRC') != -1:
                idx_test = 45
            else:
                if _dir.find('EM') != -1:
                    idx_test = 7 # 2d
                    #idx_test = 11 # 2d
                    if len(lines) > 0:
                        line_split = lines[1].split()
                        idx_test = len(line_split) - 1
                else:
                    idx_test = 9 # 3d
        else:
            idx_test = sorting_index

        scores = []
        iterationId = 0 # starts at 0 to skip first line that contains header information
        iterationId_max_score = 0
        for l in lines:
            line_split = l.split()
            if len(line_split) >= idx_test and isFloat(line_split[idx_test]):
                score = float(line_split[idx_test])
                scores.append(score)
                idx_score = idx_score + 1
                if score > max_score:
                    max_score = score
                    filename_max_score = filename
                    iterationId_max_score = iterationId
            iterationId = iterationId + 1
            if max_iteration_to_show != -1 and (iterationId*stepForOutputFiles) > max_iteration_to_show:
                break
        f.close()

        if round_max_score:
            print 'max_score=' + str(max_score)
            max_score = trunc(max_score, 2)
            print 'max_score=' + str(max_score)
            # re-iterate through the list to pick first index corresponding to the max value.
            for iterationId in range(0, len(scores)):
                if scores[iterationId] > max_score:
                    iterationId_max_score = iterationId
                    break
                
        label_name = _dir
        labels.append(label_name)
        list_max_scores.append(max_score)
        list_step_ids.append(iterationId_max_score)
        list_iteration_ids.append(iterationId_max_score*stepForOutputFiles)
        list_n_iterations.append(len(lines)*stepForOutputFiles)
        fig_idx = fig_idx + 1

        if showRunningTime:
            param_dir = _dir + '/parameter_vector0/'
            param_0 = param_dir + 'iteration_0.txt'
            n_param_files = '%d' % (len(os.listdir(param_dir)) - 1)
            param_last = param_dir + 'iteration_' + n_param_files + '.txt'
            dt = os.path.getmtime(param_last) - os.path.getmtime(param_0)
            dt_str = '%d' % dt
            list_times.append(dt_str)

        #print _dir + ' : ' + filename_max_score
        #print 'max_score = ' + str(max_score)
        #print max_score

#common_label = os.path.commonprefix(labels)
#for i in range(0,len(labels)):
#    labels[i] = ''.join(labels[i].split(common_label))

order = [i[0] for i in sorted(enumerate(list_max_scores), key=lambda x:x[1])]

for i in range(0,len(order)):
    idx = order[i]
    sc = '%.3f' % list_max_scores[idx]
    iterationId = '%d' % list_iteration_ids[idx]
    n_iterations = '%d' % list_n_iterations[idx]
    if showRunningTime:
        dt = list_times[idx]
    else:
        dt = ''
    print labels[idx] + '\t\t' + iterationId + '/' + n_iterations + '\t\t' + sc + '\t\t' + dt


if score_type == 2:
    #print '******************** TEST SET ********************'
    list_max_test_scores = []
    for i in range(0,len(order)):
        #idx = order[i]
        idx = i
        #sc = '%.3f' % list_max_scores[idx]
        iterationId = '%d' % list_iteration_ids[idx]

        _dir = labels[idx]

        if sorting_index == -1:
            if _dir.find('MSRC') != -1:
                idx_test = 45
            else:
                if _dir.find('EM') != -1:
                    idx_test = 7 # 2d
                    #idx_test = 11 # 2d
                    if len(lines) > 0:
                        line_split = lines[1].split()
                        idx_test = len(line_split) - 1
                else:
                    idx_test = 9 # 3d
        else:
            idx_test = sorting_index

        filename = labels[idx] + '/scores0/test_score.txt'
        if os.access(filename, os.F_OK):
            f = open(filename)
            lines = f.readlines()            
            if list_step_ids[idx] >= len(lines):
                # usually happens when algorithm is still writing results to the test file
                list_max_test_scores.append(0)
            else:
                l = lines[list_step_ids[idx]]
                tokens = l.split()
                if len(tokens) >= idx_test and isFloat(tokens[idx_test]):
                    list_max_test_scores.append(float(tokens[idx_test]))
                    #score = '%.3f' % float(tokens[idx_test])
                    #print labels[idx] + '\t\t' + iterationId + '/' + n_iterations + '\t\t' + score + '\t\t'

    # Display according to sorted order
    print '******************** TEST SET (SORTED) ********************'
    order_test = [i[0] for i in sorted(enumerate(list_max_test_scores), key=lambda x:x[1])]
    for i in range(0,len(order_test)):
        idx = order_test[i]
        sc = '%.3f' % list_max_test_scores[idx]
        iterationId = '%d' % list_iteration_ids[idx]
        n_iterations = '%d' % list_n_iterations[idx]
        if showRunningTime:
            dt = list_times[idx]
        else:
            dt = ''
        print labels[idx] + '\t\t' + iterationId + '/' + n_iterations + '\t\t' + sc + '\t\t' + dt
