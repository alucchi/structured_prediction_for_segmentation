import commands
import os
import sys

src_dir = '/cvlabdata1/home/biomed/EM/datasets/2012_07_02_rodentVolume/original/scans_008/testSet/'
src_extension = '.png'
dir_labels = 'labels/' # directory to output label files
dir_superpixels = 'slic_images/'
step_size = 10
cubeness = 14

exec_cmd = False
if len(sys.argv)>1:
    if sys.argv[1] == 'True' or sys.argv[1] == '1':
        print 'exec_cmd is set to TRUE'
        exec_cmd = True
    elif sys.argv[1] == 'False' or sys.argv[1] == '0':
        print 'exec_cmd is set to FALSE'
        #exec_cmd = False
    else:
        print 'Unknown first parameter'

# Create directory if it does not exist. 
if not os.access(dir_labels, os.F_OK):
    print 'Creating model ' + dir_labels
    if exec_cmd:
        os.mkdir(dir_labels)

# Create directory if it does not exist. 
if not os.access(dir_superpixels, os.F_OK):
    print 'Creating model ' + dir_superpixels
    if exec_cmd:
        os.mkdir(dir_superpixels)

# List files
ls_files = commands.getoutput('ls ' + src_dir + '*' + src_extension)
files = ls_files.split()

#print ls_files
#print files

for i in range(0,len(files)):
    #print files[i]

    basename, extension = os.path.splitext(files[i])
    seg_label = dir_labels + os.path.basename(basename) + str(step_size) + '_' + str(cubeness) + '.dat'

    if not os.access(seg_label, os.F_OK):
        cmd = './superpixel_test ' + files[i] + ' ' + str(step_size) + ' ' + str(cubeness)
        print cmd
        if exec_cmd:
            os.system(cmd)
    
        cmd = 'mv ' + os.path.basename(basename) + '.dat ' + seg_label
        print cmd
        if exec_cmd:
            os.system(cmd)

        cmd = 'mv ' + os.path.basename(basename) + '_slic_' + str(step_size) + '_' + str(cubeness) + '.png ' + dir_superpixels
        print cmd
        if exec_cmd:
            os.system(cmd)
