
import os

def read(filename, name):

    value = ''

    # open files
    fr = open(filename, 'r')
    n_lines = len(fr.readlines())

    # loop through each line in file
    fr.seek(0)
    for line in fr.readlines():
        #print line
        array = line.split('=')
        if array[0]==name:
            value = array[1].strip('\n')

    fr.close()
    return value

def write(filename, name, value):

    # open files
    fr = open(filename, 'r')
    filenamew = filename + 'w'
    fw = open(filenamew, 'w')
    n_lines = len(fr.readlines())

    # loop through each line in file
    fr.seek(0)

    updated = False
    for line in fr.readlines():
        #print line
        array = line.split('=')
        if array[0] == name:
            fw.write(name + '=' + str(value) + '\n')
            updated = True
        else:
            fw.write(line)

    if not updated:
        fw.write(name + '=' + str(value) + '\n')

    fr.close()
    fw.close()

    # make a backup and replace file
    filenamebak = filename + '.bak'
    cmd = 'mv ' + filename + ' ' + filenamebak
    os.system(cmd)
    cmd = 'mv ' + filenamew + ' ' + filename
    os.system(cmd)
