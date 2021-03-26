// ------------------------------------------------------------------
//   arch.h
//   Copyright (C) 2019-2021 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef ARCH_INCLUDED
#define ARCH_INCLUDED

extern void arch_initialize (const char *argv0);
extern unsigned arch_get_num_cores (void);
extern const char *arch_get_endianity (void);
extern const char *arch_get_ip_addr (const char *reason);
extern const char *arch_get_os (void);
extern bool arch_am_i_in_docker (void);
extern const char *arch_get_distribution (void);

#endif
