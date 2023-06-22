// ------------------------------------------------------------------
//   buf_list.h
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "genozip.h"

typedef struct {
    BufferP buf;  // address of Buffer structure (not the data!). Buffer is removed if removed bit is set with BL_SET_REMOVED 
    rom func;     // function in which this Buffer was first allocated (at which time it was added to the buffer_list)
    rom name;     // buffer name
    uint32_t code_line;
} BufListEnt;

// buflist modifications
extern void buflist_add_buf (VBlockP vb, BufferP buf, FUNCLINE);
extern void buflist_remove_buf (BufferP buf, FUNCLINE);
extern void buflist_move_buf (VBlockP vb, BufferP new_buf, rom new_name, ConstBufferP old_buf, FUNCLINE);
extern void buflist_update_vb_addr_change (VBlockP new_vb, ConstVBlockP old_vb);
extern void buflist_compact (VBlockP vb);
extern void buflist_sort (VBlockP vb, bool already_locked); // ahead of destroying

// iterators
extern BufListEnt *buflist_find_buf (VBlockP vb, ConstBufferP buf, FailType soft_fail);
extern void buflist_display (VBlockP vb);
extern void buflist_who_is_locked (void);

// freers & destroyers
extern void buflist_free_vb (VBlockP vb);
extern void buflist_free_ctx (VBlockP vb, ContextP ctx);
extern void buflist_destroy_vb_bufs (VBlockP vb, bool only_if_unused);
extern void buflist_destroy_private_and_context_vb_bufs (VBlockP vb);
extern void buflist_destroy_file_bufs (FileP file);
extern void buflist_destroy_bufs_by_name (rom name, bool compact_after_destroying);

// verifiers
extern bool buflist_test_overflows (VBlockP vb, rom msg);
extern void buflist_test_overflows_all_other_vb (VBlockP caller_vb, rom msg, bool force);
extern void buflist_test_overflows_all_vbs (rom msg);
extern bool buflist_locate (ConstBufferP buf, rom prefix);

// memory usage reporting
extern uint64_t buflist_get_memory_usage (void);
extern void buflist_show_memory (bool memory_full, unsigned max_threads, unsigned used_threads);
