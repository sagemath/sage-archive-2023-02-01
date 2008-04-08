import os, signal, sys, time, thread, threading, tempfile, hashlib, pickle
verbose = 10

def get_file_hash(fpath):
    f = file(fpath, mode="r")
    fcont = f.read()
    return hashlib.sha1(fcont).hexdigest()


def get_files_modified(flist):
    fdict = { }
    for filenm in flist:
        statinfo = os.stat(filenm)
        fdict[filenm] = statinfo.st_mtime
    return fdict

def get_files_hash(flist):
    fdict = { }
    for filenm in flist:
        fhash = get_file_hash(filenm)
        fdict[filenm] = fhash
    return fdict

def pickle_dicts(moddict, hashdict, deps, filenm):
    output = open(filenm, 'wb')
    pickle.dump(moddict, output, 2)
    pickle.dump(hashdict, output, 0)
    pickle.dump(deps, output, 2)
    output.close()

def unpickle_dicts(filenm):
    pkl_file = open(filenm, 'rb')
    moddict = pickle.load(pkl_file)
    hashdict = pickle.load(pkl_file)
    deps = pickle.load(pkl_file)
    pkl_file.close()
    return [moddict, hashdict, deps]

def create_modified_list(flist, cachefile):
    try:
        unpickled = unpickle_dicts(cachefile)
        moddict=unpickled[0]
        hashdict=unpickled[1]
    except:
        return flist

    newmod = get_files_modified(flist)

    modfiles = [ ]
    if verbose>100:
        print "--------------Modlist------------"
    for x in filelist:
        mod = newmod[x]
        try:
            oldmod = moddict[x]
            if oldmod!=mod:
                newhash = get_file_hash(x)
                oldhash = hashdict[x]
                if oldhash!=newhash:
                    modfiles.append(x)
                    if verbose>100:
                        print x
        except:
            modfiles.append(x)
            if verbose>100:
                print x
    if verbose>100:
        print "--------------EndModlist------------"
    return modfiles

def commit_new_dep_pickle(modfiles, deps, moddict, hashdict, filenm):
    moddict.update(get_files_modified(modfiles))
    hashdict.update(get_files_hash(modfiles))
    pickle_dicts(moddict, hashdict, deps, filenm)

def create_dep_dict(flist, funcdict, includedict, cwd):
    dict = { }
    for x in flist:
        ext = os.path.splitext(x)[1]
        try:
            files = (funcdict[ext])(x,includedict[ext])
            newfiles = [ ]
            for r in files:
                newfilename = os.path.normpath(r)
                newfiles.append(newfilename)
            newlist = newfiles
        except:
            newlist = [x]
        dict[x] = newlist
    return dict

def recurse_files(flist, depdict):
    fset = set(flist)
    outset = set()
    for x in fset:
        outset.update(depdict[x])
    if fset==outset:
        return fset
    else:
        return recurse_files(outset, depdict).union(fset)

def get_dep_dictionary(flist, funcdict, includedict, depdict, cwd):
    if depdict == None:
        depdict = create_dep_dict(flist, funcdict, includedict, cwd)
    fdict = { }
    for x in flist:
        fdict[x] = recurse_files([x], depdict)
    return fdict

def get_inverse_dep_dictionary(dict):
    invdepdict = { }
    for (x,y) in dict.iteritems():
        for z in y:
            try:
                invdepdict[z].update([x])
            except KeyError:
                invdepdict[z] = set([x])
    return invdepdict

def gen_modified_list(flist, invdict):
    modset = set([ ])
    for x in flist:
        try:
            modset.update(invdict[x])
        except KeyError:
            modset.update([x])
    return list(modset)

def detect_modified_list(files, moddict, hashdict, invdict):
    modlist = set([ ])
    for x in files:
        try:
            if moddict[x] != os.stat(x).st_mtime:
                if hashdict[x] != get_file_hash(x):
                    modlist.update([x])
        except:
            modlist.update([x])
    for x in set(modlist):
        try:
            modlist.update(invdict[x])
        except KeyError:
            pass
    if verbose>100:
        print "--------Dumping modlist-------------"
        print modlist
        print "--------------EndModlist------------"
    return list(modlist)

def get_compile_list(files, funcdict, includedict, cachefile, cwd):
    newdict = get_dep_dictionary(files, funcdict , includedict, None, cwd)
    if verbose>100:
        for (x,y) in newdict.iteritems():
            print "-------------Dependencies for-------------"
            print x
            print "------------------are---------------------"
            print y
            print "------------------------------------------"
    invdict = get_inverse_dep_dictionary(newdict)
    if verbose>100:
        for (x,y) in invdict.iteritems():
            print "------Inverse Dependencies for------------"
            print x
            print "------------------are---------------------"
            print y
            print "------------------------------------------"
    try:
        (moddict, hashdict, olddeps) = unpickle_dicts(cachefile)
    except:
        moddict = { }
        hashdict = { }
        olddeps = { }
    modfiles = detect_modified_list(files,moddict,hashdict,invdict)

    return (modfiles, newdict, invdict, moddict, hashdict, olddeps)


def write_cache(modfiles, newdict, invdict, moddict, hashdict, olddeps, cachefile):
    commit_new_dep_pickle(modfiles, newdict, moddict, hashdict, cachefile)

