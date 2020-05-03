// ------------------------------------------------------------------
//   arch.h
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef ARCH_INCLUDED
#define ARCH_INCLUDED

extern void arch_verify (void);

extern unsigned arch_get_num_cores (void);
extern const char *arch_get_ip_addr (void);
extern const char *arch_get_os (void);
extern bool arch_am_i_in_docker (void);

#endif
