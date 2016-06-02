import os.path
import sys
import time

def open_file(filename):
    return open(filename)

def write_file(buf, filename):
    fp = sys.stdout if filename == "stdout" else open(filename, 'w')
    fp.write(buf)
    fp.close()

def message(msg):
    print '[%s] %s'%(time.asctime(), msg)
