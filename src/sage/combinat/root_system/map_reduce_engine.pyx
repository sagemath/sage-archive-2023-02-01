############################
### Parlallel processing ###
############################
from functools import partial
cimport cython

#Rudimentary framework for performing MapReduce style computations
#DISCLAIMER: does not
cdef class MapReduceEngine():
    def __init__(self, mp_thresh=2500):
      self.worker_results = list()
      self.input_iter = None
      #Multiprocess size parameter
      self.mp_thresh = mp_thresh

    ################
    ### Reducers ###
    ################

    #Helper function for returning processed results back to parent process
    #Trivial reducer: simply collects objects with the same key in the worker
    def reduce_single_proc(self,proc):
        #Discard the zero polynomial
        reduced = set(self.worker_results)-set([tuple()])
        self.worker_results = list()
        return list(reduced)

    #All-to-one communication step
    def reduce_multi_process(self,reducer,worker_pool):
        collected_eqns = set()
        for child_eqns in worker_pool.imap_unordered(reducer,range(worker_pool._processes)):
            collected_eqns.update(child_eqns)
        return list(collected_eqns)

    #################
    ### MapReduce ###
    #################

    #Map fn across input_iter and reduce output to polynomial
    #If worker_pool is not provided, function maps and reduces on a single process.
    #If worker_pool is provided, the function attempts to determine whether it should
    #use multiprocessing based on the length of the input iterable. If it can't determine
    #the length of the input iterable then it uses multiprocessing with the default chunksize of 1
    #if chunksize is not explicitly provided.
    def map_triv_reduce(self,mapper,input_iter,reducer=None,worker_pool=None,chunksize=None,mp_thresh=None):
        #Compute multiprocessing parameters
        if worker_pool is not None:
            try:
                n = len(input_iter)
            except:
                n = mp_thresh + 1
            if chunksize is None:
                chunksize = n // (worker_pool._processes**2) + 1
        if mp_thresh is None:
          mp_thresh = self.mp_thresh
        no_mp = worker_pool is None or n < mp_thresh
        #Map phase. Casting Async Object blocks execution... Each process holds results
        #in its copy of fmats.temp_eqns
        if no_mp:
            list(map(mapper,input_iter))
        else:
            list(worker_pool.imap_unordered(mapper,input_iter,chunksize=chunksize))

        #Early termination in case no reducer is provided
        if reducer is None:
          return

        #Reduce phase
        if no_mp:
            results = reducer(0)
        else:
            results = self.reduce_multi_process(reducer,worker_pool)
        return results

    #Map the given mapper across the input iterable creted using
    #input_iter_setter, and then reduce using the provided reducer.
    #The input iterable is created in each worker process with minimal communication:
    #the parent only tells woker how to create iterable, but does NOT send an
    #iterator object or its contents

    ###INPUT:
    #input_iter_setter is a function that takes as a single argument an integer
    #0 <= proc <= worker_pool._processes. If worker_pool is None, then
    #input_iter_setter(0) is called
    #The function input_iter_setter should set mr_eng.input_iter to a desired value
    #mapper is a function that calls mr_eng.map_caller with appropriate arguments,
    #that are bound to the function inside each worker process
    #mapper does NOT return; it appends results to mr_eng.worker_results
    #reducer is a function that instructs the engine how to collect results
    #corresponding to the same key in each worker process
    #input_iter_setter, mapper, AND reducer are instructions that will be
    #processed inside each worker in the provided pool

    ###NOTES: set_input_iter MUST be called BEFORE
    def emap_no_comm(self,input_iter_setter,mapper,reducer,worker_pool=None):
        n_proc = worker_pool._processes if worker_pool is not None else 1
        params = [(child_id, n_proc) for child_id in range(n_proc)]
        self.map_triv_reduce(input_iter_setter,params,worker_pool=worker_pool,chunksize=1,mp_thresh=0)
        return self.map_triv_reduce(mapper,params,reducer,worker_pool=worker_pool,chunksize=1,mp_thresh=0)

    @cython.wraparound(False)
    @cython.nonecheck(False)
    @cython.cdivision(True)
    cpdef map_caller(self,mp_params,mapper,extra_args=None):
        cdef int child_id, i, n_proc
        child_id, n_proc = mp_params
        #Plug in additional arguments, if any
        if extra_args is not None:
            mapper = partial(mapper,*extra_args)
        if self.input_iter is None:
            raise ValueError("set_input_iter must be used to create input iterable in each worker process")
        for i, input in enumerate(self.input_iter):
            if i % n_proc == child_id:
                res = mapper(input)
                if res:
                    self.worker_results.append(res)
        #Re-set iterator
        self.input_iter = None
