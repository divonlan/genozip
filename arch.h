// ------------------------------------------------------------------
//   arch.h
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifndef ARCH_INCLUDED
#define ARCH_INCLUDED

#ifndef DISTRIBUTION
    #define DISTRIBUTION "unknown"
#endif    

extern void arch_initialize (const char *argv0);
extern unsigned arch_get_num_cores (void);
extern const char *arch_get_endianity (void);

#define ARCH_IP_LEN 16
extern const char *arch_get_ip_addr (const char *reason);
extern const char *arch_get_user_host (void);
extern const char *arch_get_os (void);

#endif
