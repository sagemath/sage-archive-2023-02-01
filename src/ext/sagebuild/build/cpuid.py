"""     asm("cpuid
         : "=a" (words[WORD_EAX]),
           "=b" (words[WORD_EBX]),
           "=c" (words[WORD_ECX]),
           "=d" (words[WORD_EDX])
         : "a" (reg));"""