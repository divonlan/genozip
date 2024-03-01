// ------------------------------------------------------------------
//   stats.h
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"
#include "buffer.h"

typedef struct {
    Did my_did_i, st_did_i;
    bool is_dict_alias;
    int64_t txt_len, z_size;
    char name[100];
    rom type;
    StrText did_i, words, dict_words, local_words, failed_ston_words, hash, uncomp_dict, comp_dict, comp_b250, uncomp_local, comp_local;
    float pc_of_txt, pc_of_z, pc_hash_occupancy, pc_ston_hash_occup;
    rom bcodec, lcodec, dcodec;
} StatsByLine;

extern void stats_generate (void);          
extern void stats_read_and_display (void);  
extern void stats_add_txt_name (rom fn);
extern void stats_finalize (void);      
    
extern void stats_add_one_program (STRp(prog_name));
extern void stats_remove_data_after_program_name (rom program);
extern rom stats_find_in_programs (rom signature);
static inline bool stats_is_in_programs (rom signature) { return (bool)stats_find_in_programs (signature); }

extern Buffer stats_programs;
