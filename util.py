import os.path
import time

def open_file(filename):
    return open(filename)

def write_file(buf, filename):
    fp = open(filename, 'w')
    fp.write(buf)
    fp.close()

def message(msg):
    print '[%s] %s'%(time.asctime(), msg)
