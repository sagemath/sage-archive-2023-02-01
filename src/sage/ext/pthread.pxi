cdef extern from "pthread.h":
    ctypedef int pthread_t       # actual type will be used by Pyrex
    ctypedef int pthread_attr_t
    int pthread_create(pthread_t *thread,
                       pthread_attr_t *attr,
                       void *(*start_routine)(void *),
                       void *arg)
    int pthread_join(pthread_t thread, void **value_ptr)
    void pthread_exit(void *value_pt)
    int pthread_attr_init(pthread_attr_t *)
    int pthread_attr_destroy(pthread_attr_t *)
    int pthread_attr_setdetachstate(pthread_attr_t *,int)
    int pthread_join(pthread_t ,void **)

    cdef enum:
       PTHREAD_CREATE_JOINABLE
