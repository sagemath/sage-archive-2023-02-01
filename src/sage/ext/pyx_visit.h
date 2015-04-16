/* Cython-compatible version of Py_VISIT */

#define Pyx_VISIT(op)                                                   \
    do {                                                                \
        if (op) {                                                       \
            int vret = __pyx_v_visit((PyObject *)(op), __pyx_v_arg);    \
            if (vret)                                                   \
                return vret;                                            \
        }                                                               \
    } while (0)

