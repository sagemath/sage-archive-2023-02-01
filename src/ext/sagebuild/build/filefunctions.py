import os, signal, sys, time, thread, threading, tempfile

def safemkdirs(path):
    try:
        os.makedirs(path)
    except OSError:
        pass

def safesymlink(target, path):
    try:
        os.symlink(target,path)
    except OSError:
        pass