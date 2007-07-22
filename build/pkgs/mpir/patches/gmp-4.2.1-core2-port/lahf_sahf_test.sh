cat > tmp_lahf_sahf_test.c <<EOF
#include<stdio.h>
int are_available_lahf_sahf(void)
{
register unsigned long long rcx asm ("rcx");
rcx = 0;
asm ("movl\t\$0x80000001,%eax");
asm ("cpuid");
return(rcx & 0x0000000000000001LL);
}

main()
{
if (are_available_lahf_sahf())
{
printf("Yes");
}
else
{
printf("No");
}
}
EOF
gcc -m64 tmp_lahf_sahf_test.c -o tmp_lahf_sahf_test > /dev/null 2>&1
IS_LAHF_AVAIL=`./tmp_lahf_sahf_test`;
if [ $IS_LAHF_AVAIL == "Yes" ] ; then
    echo -n "Yes";
else
    echo -n "No";
fi
rm -f tmp_lahf_sahf_test.c
rm -f tmp_lahf_sahf_test

