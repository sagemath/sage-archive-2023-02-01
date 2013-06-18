
#if defined(__APPLE__)

// Returns the virtual size of the current process. This will match what "top" reports for VSIZE
unsigned long long darwin_virtual_size();

#endif

