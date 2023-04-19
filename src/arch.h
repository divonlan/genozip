// ------------------------------------------------------------------
//   arch.h
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

extern void arch_initialize (rom argv0);
extern unsigned arch_get_num_cores (void);
extern double arch_get_physical_mem_size (void);
extern rom arch_get_endianity (void);
extern void arch_set_locale (void);

#define ARCH_IP_LEN 16
extern rom arch_get_ip_addr (rom reason);
extern rom arch_get_host (void);
extern rom arch_get_user_host (void);
extern rom arch_get_os (void);
extern rom arch_get_distribution (void);
extern rom arch_get_executable (FailType soft_fail);
extern rom arch_get_argv0 (void);
extern bool arch_is_valgrind (void);

extern Timestamp arch_timestamp (void);

static inline uint32_t arch_time_lap (uint128_t ts_start) // in msec
{ return (uint32_t)((arch_timestamp() - ts_start) / 1000000); }
