import os

_proc_status = '/proc/%d/status' % os.getpid()
#_scale = {'kB': 1024.0, 'mB': 1024.0*1024.0,
#          'KB': 1024.0, 'MB': 1024.0*1024.0}

_scale = {'kB': 1, 'mB': 1024, 'gB': 1024 * 1024,
          'KB': 1, 'MB': 1024, 'GB': 1024 * 1024}


def _VmB(VmKey):
    '''Private.
    '''
    global _proc_status, _scale
     # get pseudo file  /proc/<pid>/status
    try:
        t = open(_proc_status)
        v = t.read()
        t.close()
    except:
        return float('nan')  # non-Linux?
     # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
    i = v.index(VmKey)
    v = v[i:].split(None, 3)  # whitespace
    if len(v) < 3:
        return float('nan')  # invalid format?
     # convert Vm value to bytes
  #  return float(v[1]) * _scale[v[2]]
    return int(v[1]) * _scale[v[2]]


def memory(since=0):
    '''Return memory usage in kilobytes.
    '''
    return _VmB('VmSize:') - since


def resident(since=0):
    '''Return resident memory usage in kilobytes.
    '''
    return _VmB('VmRSS:') - since


def memorypeak(since=0):
    '''Return memory usage peak in kilobytes.
    '''
    try:
        return _VmB('VmPeak:') - since
    except:
        return float('nan')  # old Linux?


def residentpeak(since=0):
    '''Return resident memory usage peak in kilobytes.
    '''
    try:
        return _VmB('VmHWM:') - since
    except:
        return float('nan')  # old Linux?


def stacksize(since=0):
    '''Return stack size in kilobytes.
    '''
    return _VmB('VmStk:') - since
