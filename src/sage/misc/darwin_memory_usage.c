
#if defined(__APPLE__)

#include <mach/host_info.h>
#include <mach/mach_host.h>
#include <mach/shared_region.h>
#include <mach/mach.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/sysctl.h>

static inline cpu_type_t cpu_type() {
    cpu_type_t cputype = 0;

    int mib[CTL_MAXNAME];
    size_t len = CTL_MAXNAME;
    if (sysctlnametomib("sysctl.proc_cputype", mib, &len) != -1) {
        mib[len] = getpid();
        len++;

        size_t cputypelen = sizeof(cputype);
        if (sysctl(mib, len, &cputype, &cputypelen, 0, 0) == -1) {
            cputype = 0;
        }
    }

    return cputype;
}

static inline mach_vm_size_t shared_region_size(cpu_type_t cputype) {
    switch(cputype) {
        case CPU_TYPE_POWERPC:
            return SHARED_REGION_SIZE_PPC;
        case CPU_TYPE_POWERPC64:
            return SHARED_REGION_SIZE_PPC64;
        case CPU_TYPE_I386:
            return SHARED_REGION_SIZE_I386;
        case CPU_TYPE_X86_64:
            return SHARED_REGION_SIZE_X86_64;
        default: // unknown CPU type
            return 0;
    }
}

static inline int in_shared_region(cpu_type_t cputype, mach_vm_address_t address) {
    if (cputype == CPU_TYPE_POWERPC &&
        address >= SHARED_REGION_BASE_PPC &&
        address <= (SHARED_REGION_BASE_PPC + SHARED_REGION_SIZE_PPC)) {
        return 1;
    }

    if (cputype == CPU_TYPE_POWERPC64 &&
        address >= SHARED_REGION_BASE_PPC64 &&
        address <= (SHARED_REGION_BASE_PPC64 + SHARED_REGION_SIZE_PPC64)) {
        return 1;
    }

    if (cputype == CPU_TYPE_I386 &&
        address >= SHARED_REGION_BASE_I386 &&
        address <= (SHARED_REGION_BASE_I386 + SHARED_REGION_SIZE_I386)) {
        return 1;
    }

    if (cputype == CPU_TYPE_X86_64 &&
        address >= SHARED_REGION_BASE_X86_64 &&
        address <= (SHARED_REGION_BASE_X86_64 + SHARED_REGION_SIZE_X86_64)) {
        return 1;
    }

    return 0;
}

unsigned long long darwin_virtual_size()
{
    kern_return_t error;
    task_t task;
    struct task_basic_info_64 taskinfo;
    cpu_type_t cputype;
    mach_msg_type_number_t count;
    mach_vm_size_t size;
    mach_vm_address_t address;
    mach_port_t object_name;
    vm_region_top_info_data_t info;
    mach_vm_size_t	vsize;
    mach_vm_size_t	empty;
    int has_shared_regions;

    empty = 0;

    count = TASK_BASIC_INFO_64_COUNT;
    task = mach_task_self();
    error = task_info(task, TASK_BASIC_INFO_64, (task_info_t)&taskinfo, &count);

    if (error != KERN_SUCCESS) {
        return 0;
    }

    vsize = taskinfo.virtual_size;

    cputype = cpu_type();

    // Go through all the vm regions and check to see if we should count them in the vsize or not
    for (address = 0, has_shared_regions = 0; ; address += size) {
        count = VM_REGION_TOP_INFO_COUNT;
        if (mach_vm_region(task, &address, &size, VM_REGION_TOP_INFO, (vm_region_info_t)&info, &count, &object_name) != KERN_SUCCESS) {
            // There are no more vm regions to look at.
            break;
        }

        if (in_shared_region(cputype, address)) {
            // Check if this process has the globally shared text and data regions mapped in.
            // If so, set has_shared_regions to 1 and so we only check once.
            if (has_shared_regions == 0 && info.share_mode == SM_EMPTY) {
                vm_region_basic_info_data_64_t basic_info;

                count = VM_REGION_BASIC_INFO_COUNT_64;
                if (mach_vm_region(task, &address, &size, VM_REGION_BASIC_INFO, (vm_region_info_t)&basic_info, &count, &object_name) != KERN_SUCCESS) {
                    break;
                }

                if (basic_info.reserved) {
                    has_shared_regions = 1;
                }
            }

            // Skip the vm region if it is not a shared private region.
            if (info.share_mode != SM_PRIVATE) {
                continue;
            }
        }

        if (info.share_mode == SM_EMPTY) {
            empty += size;
        }
    }

    // Subtract out the globally shared text and data region.
    if (has_shared_regions == 1) {
        vsize -= shared_region_size(cputype);
    }

    // Subtract out the empty pages (pagezero, stack guard, etc)
    vsize -= empty;

    return vsize;
}

#endif
