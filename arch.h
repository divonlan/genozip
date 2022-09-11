// ------------------------------------------------------------------
//   arch.h
//   Copyright (C) 2019-2022 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is propeitary, not open source software. Modifying the source code is strictly not permitted
//   and subject to penalties specified in the license.

#pragma once

extern void arch_initialize (rom argv0);
extern unsigned arch_get_num_cores (void);
extern rom arch_get_endianity (void);

#define ARCH_IP_LEN 16
extern rom arch_get_ip_addr (rom reason);
extern rom arch_get_user_host (void);
extern rom arch_get_os (void);
extern rom arch_get_distribution (void);
extern bool arch_is_wsl (void);
extern rom arch_get_executable (rom argv0);

