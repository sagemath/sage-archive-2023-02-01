// On some systems, macros "minor()" and "major()" are defined in system header files. This will undefine those:
#ifdef major
    #undef major
#endif
#ifdef minor
    #undef minor
#endif
