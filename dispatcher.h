// ------------------------------------------------------------------
//   dispatcher.h
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#ifndef DISPATCHER_INCLUDED
#define DISPATCHER_INCLUDED

#include "genozip.h"
#include "buffer.h"

typedef void *Dispatcher;
extern Dispatcher dispatcher_init (unsigned max_threads, unsigned previous_vb_i, FileP vcf_file, FileP z_file,
                                   bool test_mode, bool is_last_file, const char *filename);
extern void dispatcher_pause (Dispatcher dispatcher);
extern void dispatcher_resume (Dispatcher dispatcher, FileP vcf_file);
extern void dispatcher_finish (Dispatcher *dispatcher, unsigned *last_vb_i);

typedef void (*DispatcherFuncType)(VariantBlockP);
extern void dispatcher_compute (Dispatcher dispatcher, DispatcherFuncType func);
extern VariantBlockP dispatcher_generate_next_vb (Dispatcher dispatcher);       
extern bool dispatcher_has_processed_vb (Dispatcher dispatcher, bool *is_final);                                  
extern VariantBlockP dispatcher_get_processed_vb (Dispatcher dispatcher, bool *is_final);
extern bool dispatcher_has_free_thread (Dispatcher dispatcher);
extern VariantBlockP dispatcher_get_next_vb (Dispatcher dispatcher);
extern void dispatcher_finalize_one_vb (Dispatcher dispatcher, ConstFileP file, long long vcf_data_written_so_far, uint64_t bytes_compressed);
extern void dispatcher_input_exhausted (Dispatcher dispatcher);
extern bool dispatcher_is_done (Dispatcher dispatcher);
extern bool dispatcher_is_input_exhausted (Dispatcher dispatcher);

#endif
