// ------------------------------------------------------------------
//   fast_shared.c
//   Copyright (C) 2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include "fast_private.h"
#include "file.h"
#include "endianness.h"
#include "piz.h"
#include "reconstruct.h"
#include "dict_id.h"
#include "mutex.h"
#include "strings.h"
#include "regions.h"
#include "random_access.h"

unsigned fast_vb_size (void) { return sizeof (VBlockFAST); }
unsigned fast_vb_zip_dl_size (void) { return sizeof (ZipDataLineFAST); }

void fast_vb_release_vb (VBlockFAST *vb)
{
    vb->last_line = 0;
    vb->contig_grepped_out = false;
    vb->pair_num_lines = vb->pair_vb_i = 0;
    
    if (z_file && (z_file->data_type == DT_FASTA || z_file->data_type == DT_REF))
        vb->contexts[FASTA_NONREF].local.len = 0; // len might be is used even though buffer is not allocated (in make-ref)

    buf_free (&vb->genobwa_show_line);

    FREE (vb->optimized_desc);
    vb->optimized_desc_len = 0;
}

void fast_vb_destroy_vb (VBlockFAST *vb)
{
    buf_destroy (&vb->genobwa_show_line);
}


// called by I/O thread in fast_piz_read_one_vb, in case of --grep, to decompress and reconstruct the desc line, to 
// see if this vb is included. 
bool fast_piz_test_grep (VBlockFAST *vb)
{
    ARRAY (const unsigned, section_index, vb->z_section_headers);

    SectionHeaderVbHeader *header = (SectionHeaderVbHeader *)(vb->z_data.data + section_index[0]);
    vb->first_line       = BGEN32 (header->first_line);
    vb->lines.len        = BGEN32 (header->num_lines);
    vb->vb_data_size     = BGEN32 (header->vb_data_size);
    vb->longest_line_len = BGEN32 (header->longest_line_len);

    // in case of --unbind, the vblock_i in the 2nd+ component will be different than that assigned by the dispatcher
    // because the dispatcher is re-initialized for every sam component
    if (flag.unbind) vb->vblock_i = BGEN32 (header->h.vblock_i);
    
    // we only need room for one line for now 
    buf_alloc (vb, &vb->txt_data, vb->longest_line_len, 1.1, "txt_data");

    // uncompress & map desc field (filtered by piz_is_skip_section)
    vb->grep_stages = GS_TEST; // tell piz_is_skip_section to skip decompressing sections not needed for determining the grep
    piz_uncompress_all_ctxs ((VBlockP)vb, 0);
    vb->grep_stages = GS_UNCOMPRESS; // during uncompress in the compute thread, uncompress only what was not already uncompressed here

    // reconstruct each description line and check for string matching with flag.grep
    bool found = false, match = false;

    Context *desc_ctx =  &vb->contexts[vb->data_type == DT_FASTQ ? FASTQ_DESC : FASTA_DESC];
    desc_ctx->iterator.next_b250 = FIRSTENT (uint8_t, desc_ctx->b250); 

    uint32_t num_descs; // number of desciption lines in this VB
    if (vb->data_type == DT_FASTQ) {
        vb->line_i = 4 * vb->first_line;
        num_descs = (uint32_t)vb->lines.len; // every read has a description line
    }
    else {
        vb->line_i = vb->first_line;
        num_descs = random_access_num_chroms_start_in_this_vb (vb->vblock_i);
    }

    // iterate on all DESCs of this VB
    for (uint32_t desc_i=0; desc_i < num_descs; desc_i++) { 

        reconstruct_from_ctx (vb, desc_ctx->did_i, 0, true);

        *AFTERENT (char, vb->txt_data) = 0; // terminate the desc string

        match = flag.grep && !!strstr (vb->txt_data.data, flag.grep); // note: this function is also called due to --regions in FASTA

        if (!match && flag.regions && vb->data_type == DT_FASTA) 
            match = fasta_piz_is_grepped_out_due_to_regions (vb, vb->txt_data.data);

        vb->txt_data.len = 0; // reset

        if (match) { // 
            found = true; // we've found a match to the grepped string
            if (vb->data_type == DT_FASTQ) break; // for FASTA, we need to go until the last line, for FASTQ, we can break here
        }

        if (vb->data_type == DT_FASTQ) vb->line_i += 4; // note: for FASTA we have no idea what txt line we're on, because we're only tracking DESC lines
    }

    // last FASTA - carry over whether its grepped to the next VB - in case next VB starts not from the description line
    // similarly, note whether the previous VB ended with a grepped sequence. If previous VB didn't have any description
    // i.e the entire VB was a sequence that started in an earlier VB - the grep status of the easier VB is carried forward
    if (vb->data_type == DT_FASTA) 
        // if the last contig of the previous vb was grepped in - then include this VB anyway
        found = fasta_piz_initialize_contig_grepped_out (vb, desc_ctx->b250.len > 0, match) || found;

    // reset iterators - piz_fast*_reconstruct_vb will use them again 
    ctx_init_iterator (desc_ctx);
    for (DidIType did_i=0; did_i < vb->num_contexts; did_i++) {
        
        Context *ctx = &vb->contexts[did_i];
        if (dict_id_is_type_1 (ctx->dict_id)) {
            ctx_init_iterator (ctx);
            ctx->last_delta = ctx->last_value.f = 0;
        }
    }

    return found; 
}

