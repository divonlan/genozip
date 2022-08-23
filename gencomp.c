// ------------------------------------------------------------------
//   gencomp.c - "generated component"
//   Copyright (C) 2021-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

#include <errno.h>
#include <pthread.h>
#include "libdeflate/libdeflate.h"
#include "genozip.h"
#include "file.h"
#include "gencomp.h"
#include "zip.h"
#include "sections.h"
#include "sam.h"
#include "fastq.h"
#include "codec.h"
#include "bgzf.h"

//-----------------------
// Types & macros
//-----------------------

#define MAX_QUEUE_SIZE_BITS 16 // fields in queue and QBits are 16 bits
#define MAX_QUEUE_SIZE ((1 << MAX_QUEUE_SIZE_BITS)-2)
#define END_OF_LIST (MAX_QUEUE_SIZE+1)

typedef struct {
    BufferP gc_txts;         // Buffer array - each buffer is either on the "ready-to-dispatch queue" or on the "unused buffers' stack"
    uint16_t queue_size;     // number of buffers allocated on queue
    volatile uint16_t queue_len;  // number of buffers ready for dispatching (NOT including buffers offloaded to disk)
    uint16_t next_to_leave;  // Head of queue (FIFO): doubly-linked list of buffers 
    uint16_t last_added;     // Tail of queue 
    uint16_t next_unused;    // Top of stack: Unused gc_txts buffers stack (LIFO): linked-list of buffers
} QueueStruct;

// overlays Buffer.param of queue[gct].gc_txts[buf_i]
typedef struct __attribute__ ((__packed__)) {
    uint32_t num_lines : 27; // VB size is limited to 2GB and VCF spec mandates at least 16 characters per line (8 mandatory fields)
    uint32_t comp_i    : 2;
    uint32_t unused    : 3;
    uint16_t next;           // queue: towards head ; stack: down stack
    uint16_t prev;           // queue: towards tail ; stack: always END_OF_LIST
} QBits;

typedef union {
    QBits bits;
    int64_t value;
} QBitsType;
#define GetQBit(buf_i,field) ((QBitsType){ .value = queueP[gct].gc_txts[buf_i].param }.bits.field)
#define SetQBit(buf_i,field) ((QBitsType*)&queueP[gct].gc_txts[buf_i].param)->bits.field
#define SetQBits(buf_i)      ((QBitsType*)&queueP[gct].gc_txts[buf_i].param)->bits

typedef struct {
    uint32_t comp_len;
    uint32_t num_lines;
    CompIType comp_i;
} TxtDataInfoType;

typedef struct {
    char *name;          
    FILE *fp;   
    Buffer offload_info; // array of TxtDataInfoType - one entry per txt_data buffer written to the file. 
    Buffer thread_data;
    Buffer thread_data_comp;
    bool has_thread;
    pthread_t thread;
} DepnStruct;

typedef struct  {
    GencompType type;
    Buffer txt_data;     // we absorb gencomp lines separately for each component, and then flush a single-component VB to the queue - the queue may have VBs from several components
    int64_t num_lines;   // for stats
} CompStruct;

// Thread safety model:
// At most 2 threads can access our data structures concurrently:
// "Absorb thread" - this is either one of the compute threads (selected using a serializing VB) or the main thread 
// for the final flush, after all compute threads are done absorbing. This thread accesses only the "Absorbing functions".
// "Dispatcher thread" - called by the main thread from zip_prepare_one_vb_for_dispatching, calling "Dispatcher functions".
// After all VBs are absorbed, a final call to flush occurs from the main thread, it is the only thread running and is considered
// to be both "Absorb thread" and "Dispatcher thread".
//
// gc_protected protects variables that may be accessed by both "Absorbing functions" and "Dispatcher functions" running concurrently.
static Mutex gc_protected = {};
static Serializer gc_vb_serializer = {};

// --------------------------------------------------------------------------------------
// unprotected data structures - accessed only by "Dispatcher functions" (main thread)
// --------------------------------------------------------------------------------------
static bool sam_finished_ingesting_prim = false;
static VBIType num_SAM_PRIM_vbs_ingested = 0;
static VBIType num_vbs_dispatched[NUM_GC_TYPES] = {}; // by gencomp type
static DepnStruct depn = {};              // doesn't need protection because accessed by Dispatcher functions only after Absorption is complete

static enum { DEPN_NONE, DEPN_OFFLOAD, DEPN_REREAD } depn_method; // constant after initialization in gencomp_initialize

// --------------------------------------------------------------------------------------
// Mutex-protected data structures - accessed by both "Dispatcher functions" (main thread)
// and "Absorbing functions" (serialized compute thread) running concurrently (P final)
// --------------------------------------------------------------------------------------
static bool finished_absorbingP = false;
static VBIType num_MAIN_vbs_absorbedP = 0; 
static QueueStruct queueP[NUM_GC_TYPES] = {}; // queue of txt_data's [1] out-of-band (used for SAM PRIM, DVCF PRIM + LUFT) [2] DEPN (used for SAM DEPN)
static CompStruct componentsP[MAX_GEN_COMP+1] = {};

// --------------------------------------------------------------------------------------
// unprotected data is accessed by absorbing threads until absorbing is done, 
// and after by dispatcher functions
// --------------------------------------------------------------------------------------
static Buffer reread_depn_lines = {}; // array of type RereadLine

//--------------------------------------------------
// Seg: adding gencomp lines to vb->gencomp
//--------------------------------------------------

// ZIP compute thread: store location where this gc line should be inserted    
void gencomp_seg_add_line (VBlockP vb, CompIType comp_i, STRp(line)/*pointer into txt_data*/)
{
    ASSERT (line_len <= GENCOMP_MAX_LINE_LEN, "line_len=%u is beyond maximum of %u", line_len, GENCOMP_MAX_LINE_LEN);
    
    buf_alloc (VB, &vb->gencomp_lines, 1, 5000, GencompLineIEntry, 2, "gencomp_lines");
 
    BNXT (GencompLineIEntry, vb->gencomp_lines) = 
        (GencompLineIEntry){ .line_i     = vb->line_i, 
                             .line_index = BNUMtxt (line),
                             .line_len   = line_len,
                             .comp_i     = comp_i };

    // If we're might re-read depn lines from the txt file, we store their coordinates in the txt file
    if (componentsP[comp_i].type == GCT_DEPN && depn_method == DEPN_REREAD) {
        if (txt_file->codec == CODEC_BGZF) {
            uint64_t bb_i = vb->vb_bgzf_i + vb->bgzf_blocks.current_bb_i;
            ASSERT (bb_i <= MAX_BB_I, "BGZF bb_i=%"PRIu64" exceeds maximum of %"PRIu64, bb_i, MAX_BB_I);

            BLST (GencompLineIEntry, vb->gencomp_lines)->offset = (LineOffset){ .bb_i = bb_i, .uoffset = vb->line_bgzf_uoffset };
        }
        else // CODEC_NONE
            BLST (GencompLineIEntry, vb->gencomp_lines)->offset = (LineOffset){ .offset = vb->vb_position_txt_file + BNUMtxt (line) };
    }
}

//--------------------------------------------------
// Misc functions
//--------------------------------------------------

// true if VB belongs to a generated componentsP
bool gencomp_comp_eligible_for_digest (VBlockP vb)
{
    CompIType comp_i = vb ? vb->comp_i : flag.zip_comp_i;
    
    return (comp_i == COMP_MAIN) || 
           (vb && vb->data_type == DT_FASTQ) || (!vb && z_file->data_type == DT_FASTQ); // FASTQ comp_i==1 is R2 and not a generated component
}

static void debug_gencomp (rom msg, bool needs_lock)
{
    if (needs_lock) mutex_lock (gc_protected); 

    QueueStruct *q1 = &queueP[GCT_OOB], *q2 = &queueP[GCT_DEPN];

    iprintf ("%-12.12s OOB: size=%u len=%u tail=%-2d head=%-2d unused=%-2d DEPN: S=%u L=%u T=%-2d H=%-2d U=%-2d #on_disk=%u\n",
             msg, 
             q1->queue_size, q1->queue_len, (int16_t)q1->last_added, (int16_t)q1->next_to_leave, (int16_t)q1->next_unused,
             q2->queue_size, q2->queue_len, (int16_t)q2->last_added, (int16_t)q2->next_to_leave, (int16_t)q2->next_unused, (int)depn.offload_info.len);

    if (needs_lock) mutex_unlock (gc_protected);
}

//-----------------------------------------------------------------------------
// Initalization and finalization - called by main thread outside of dispatcher
// loop - no compute threads are running
//-----------------------------------------------------------------------------

// main thread
void gencomp_initialize (CompIType comp_i, GencompType gct) 
{
    // case: first call for this z_file
    componentsP[COMP_MAIN] = (CompStruct){ .type = GCT_NONE }; 
    componentsP[comp_i]    = (CompStruct){ .type = gct      };

    if (!gc_vb_serializer.initialized) { // initialize once for all components
        serializer_initialize (gc_vb_serializer);
        mutex_initialize (gc_protected);
    }

    // add to buffer list. we can't allocate yet because segconf.vb_size is not known yet
    buf_add_to_buffer_list (evb, &componentsP[comp_i].txt_data);

    // initialize queue. note: same-type components share a queue: in DVCF, LUFT and PRIM are both out-of-band.
    #if MAX_QUEUE_SIZE < MAX_GLOBAL_MAX_THREADS 
    #error MAX_GLOBAL_MAX_THREADS must be at most MAX_QUEUE_SIZE or else the queue_size might overflow
    #endif

    if (gct == GCT_DEPN) {

        buf_add_to_buffer_list_(evb, &depn.thread_data, "depn.thread_data");
        buf_add_to_buffer_list_(evb, &depn.thread_data_comp, "depn.thread_data_comp");

        // if we cannot re-read the depn lines from the file, we will offload them to disk
        if ((txt_file->codec != CODEC_BGZF && txt_file->codec != CODEC_NONE) || txt_file->redirected || txt_file->is_remote) {
            depn_method = DEPN_OFFLOAD;

            depn.name = MALLOC (strlen (z_file->name) + 20);
            sprintf (depn.name, "%s.DEPN", z_file->name);
            buf_add_to_buffer_list (evb, &depn.offload_info);
        }

        else {
            depn_method = DEPN_REREAD;
            buf_add_to_buffer_list_(evb, &reread_depn_lines, "reread_depn_lines");
        }
    }

    if (!queueP[gct].gc_txts) {
        queueP[gct] = (QueueStruct) {
            .queue_size    = (gct == GCT_OOB)             ? global_max_threads // OOB:  global_max_threads is the theoretical maximum number of gencomp txt_data generatable from global_max_threads concurrently running ZIP compute threads
                           : (depn_method == DEPN_REREAD) ? global_max_threads // DEPN: when re-reading from txt_file data that doesn't fit into the queue - a larger queue is beneficial
                           :                                1,                 // DEPN: when writing a file with the data that doesn't fit into the queue - No strong benefit of having more buffers in memory, see: private/internal-docs/gencomp-depn-memory-queue-vs-disk.txt
            .next_unused   = 0,           // 0 means gc_txts[0]
            .last_added    = END_OF_LIST, // ready-to-dispatch queue is initially empty
            .next_to_leave = END_OF_LIST
        };

        queueP[gct].gc_txts = CALLOC (sizeof (Buffer) * queueP[gct].queue_size);

        // add to evb buf_list, so we can buf_alloc in other threads (similar to INIT in file.c)
        for (int i=0; i < queueP[gct].queue_size; i++) 
            buf_add_to_buffer_list_(evb, &queueP[gct].gc_txts[i], "queueP.gc_txts");

        // add all buffers to "unused stack"
        for (int i=0; i < queueP[gct].queue_size-1; i++) {
            SetQBit (i, next) = i+1;
            SetQBit (i, prev) = END_OF_LIST; // always EOL for stack
        }

        SetQBit (queueP[gct].queue_size-1, next) = END_OF_LIST;
    }
}

// main thread
void gencomp_destroy (void)
{
    for (GencompType gct=1; gct < NUM_GC_TYPES; gct++)
        if (queueP[gct].gc_txts) {
            for (int i=0; i < queueP[gct].queue_size; i++) 
                buf_destroy (queueP[gct].gc_txts[i]);
            FREE(queueP[gct].gc_txts);
        }

    for (CompIType comp_i=1; comp_i <= MAX_GEN_COMP; comp_i++)
        buf_destroy (componentsP[comp_i].txt_data);

    if (depn_method == DEPN_OFFLOAD) {
        if (depn.fp) {
            fclose (depn.fp);
            remove (depn.name);
        }
        FREE (depn.name);
        buf_destroy (depn.offload_info);
        buf_destroy (depn.thread_data);
        buf_destroy (depn.thread_data_comp);
    }

    else if (depn_method == DEPN_REREAD) {
        buf_destroy (reread_depn_lines);
    }

    serializer_destroy (gc_vb_serializer);
    mutex_destroy (gc_protected);
    memset (queueP,   0, sizeof (queueP));
    memset (componentsP, 0, sizeof (componentsP));
    memset (&depn, 0, sizeof (depn));
    memset ((void*)num_vbs_dispatched, 0, sizeof(num_vbs_dispatched));
    finished_absorbingP = sam_finished_ingesting_prim = false;
    num_SAM_PRIM_vbs_ingested = num_MAIN_vbs_absorbedP = 0;
    depn_method = DEPN_NONE;
}

// main thread (from stats)
int64_t gencomp_get_num_lines (CompIType comp_i)
{
    return componentsP[comp_i].num_lines;
}

//-----------------------------------------------------------------------------------------------
// "Absorbing functions" - Absorbing data from VBs into gencomp queues - these are called by
// the Absorb thread (for the final flush, the Absorbing and Dispatching is the same thread)
//-----------------------------------------------------------------------------------------------

static uint32_t compress_depn_buf (BufferP comp_buf)
{
    static VBlockP compress_depn_vb = NULL;

    if (!compress_depn_vb) 
        // get VB - needed only for its codec_bufs
        compress_depn_vb = vb_initialize_nonpool_vb (VB_ID_COMPRESS_DEPN, DT_BAM, "compress_depn_buf");
    
    uint32_t uncomp_len = depn.thread_data.len32;
    uint32_t comp_len = codec_RANB_est_size (CODEC_RANS8, uncomp_len);

    buf_alloc (evb, comp_buf, 0, comp_len + sizeof (uint32_t), char, 1, NULL); // added to evb's buf_list before
    *B1ST32 (*comp_buf) = uncomp_len;

    // about 3X compression (BAM) with no impact on time (compression+decompression time is a wash with the saving of I/O time)
    codec_RANB_compress (compress_depn_vb, NULL, NULL, depn.thread_data.data, &uncomp_len, NULL, comp_buf->data + sizeof(uint32_t), 
                         &comp_len, false, NULL);
    
    return (comp_buf->len = comp_len + sizeof (uint32_t));
}

static void *gencomp_do_offload (void *info_)
{
    TxtDataInfoType *info = (TxtDataInfoType *)info_;
    uint32_t uncomp_len = depn.thread_data.len32;
    info->comp_len = compress_depn_buf (&depn.thread_data_comp);

    ASSERT (1 == fwrite (STRb(depn.thread_data_comp), 1, depn.fp), 
            "Failed to write %"PRIu64" bytes to temporary file %s: %s", depn.thread_data_comp.len, depn.name, strerror (errno));

    if (flag.debug_gencomp) 
        iprintf ("Wrote to disk: buf=%u num_lines=%u uncomp_len=%u comp_len=%u uncomp_alder32=%u comp_adler32=%u\n",
                 BNUM (depn.offload_info, info), info->num_lines, uncomp_len, info->comp_len, 
                 adler32 (1, STRb(depn.thread_data)), adler32 (1, STRb(depn.thread_data_comp))); 

    depn.thread_data_comp.len = depn.thread_data.len = 0;

    return NULL;
}

// called by gencomp_flush with gc_protected locked. 
static void gencomp_offload_DEPN_to_disk (CompIType comp_i, bool is_final_flush)
{
    // open file upon first write 
    if (!depn.fp) {
        depn.fp = fopen (depn.name, WRITEREAD); 
        ASSERT (depn.fp, "fopen() failed to open %s: %s", depn.name, strerror (errno));

        // unlink from directory so file gets deleted on close and hopefully also prevents unnecessary flushing
        unlink (depn.name); // ignore errors (this doesn't work on NTFS)
    }

    BufferP buf = &componentsP[comp_i].txt_data;

    buf_alloc (evb, &depn.offload_info, 1, 100, TxtDataInfoType, 2, "depn.offload_info");
    BNXT (TxtDataInfoType, depn.offload_info) = (TxtDataInfoType){
        .comp_i    = comp_i,
        .num_lines = buf->count
        // .comp_len is set by gencomp_do_offload
    };
    TxtDataInfoType *info = BLST (TxtDataInfoType, depn.offload_info);

    // copy data to another buffer, as this one will be freed and recycled
    buf_copy (evb, &depn.thread_data, buf, uint8_t, 0, 0, NULL);

    int err;
    if (!is_final_flush && global_max_threads > 1) {
        // compress and write to disk in separate thread, as to not delay releasing this VB. Only one thread is run in parallel as we join
        // the previous before creating a new one.
        ASSERT (!(err = pthread_create (&depn.thread, NULL, gencomp_do_offload, info)), 
                "failed to create gencomp_do_offload thread: %s", strerror(err));
        
        depn.has_thread = true;
    }
    else 
        gencomp_do_offload(info);

    if (flag.debug_gencomp) debug_gencomp ("offloaded", false);
}

static void *gencomp_compress_depn (void *comp_buf)
{
    compress_depn_buf ((BufferP)comp_buf);
    depn.thread_data_comp.len = depn.thread_data.len = 0;

    return NULL;
}

// ZIP: called with gc_protected locked.
// moves data from componentsP.txt_data where it was absorbing lines from MAIN VBs to the GetQBit or DEPN queue
static void gencomp_flush (CompIType comp_i, bool is_final_flush) // final flush of gct, not of comp_i!
{
    if (!componentsP[comp_i].txt_data.len) return; // nothing to flush

    componentsP[comp_i].num_lines += componentsP[comp_i].txt_data.count; // for stats

    GencompType gct = componentsP[comp_i].type;

    // pop a buffer from the unused stack
    uint16_t buf_i = queueP[gct].next_unused;

    // sanity
    ASSERT ((buf_i == END_OF_LIST) == (queueP[gct].queue_len == queueP[gct].queue_size), 
            "Queue is broken: buf_i=%d END_OF_LIST=%d queueP[%u].queue_len=%u queueP[%u].queue_size=%u",
            buf_i, END_OF_LIST, gct, queueP[gct].queue_len, gct, queueP[gct].queue_size);

    // wait for previous DEPN compression / offload to finish
    if (gct == GCT_DEPN && depn.has_thread) {
        int err;
        ASSERT (!(err = pthread_join (depn.thread, NULL)), "pthread_join failed: %s", strerror (err));
        depn.has_thread = false;
    }

    // if there are no more room on the queue - offload to disk (only if GCT_DEPN)
    if (buf_i == END_OF_LIST) {
        // case: Out-of-band: this is an error
        ASSERT (gct != GCT_OOB, "Out-of-band gencomp queue is full: gc_queue_size=%d", queueP[gct].queue_size);

        // case: DEPN - we offload one buffer from the queue to the file
        gencomp_offload_DEPN_to_disk (comp_i, is_final_flush);
    }

    else {
        queueP[gct].next_unused = is_final_flush ? END_OF_LIST : GetQBit(buf_i, next);

        // add the buffer at the tail of the queue
        if (queueP[gct].last_added != END_OF_LIST)
            SetQBit (queueP[gct].last_added, prev) = buf_i;

        SetQBits(buf_i) = (QBits){ .next      = queueP[gct].last_added,
                                   .prev      = END_OF_LIST,
                                   .num_lines = componentsP[comp_i].txt_data.count,
                                   .comp_i    = comp_i } ;
        
        queueP[gct].last_added = buf_i;
        queueP[gct].queue_len++;

        if (queueP[gct].next_to_leave == END_OF_LIST) {
            ASSERT0 (queueP[gct].queue_len==1, "Expecting that if queue's head==tail, then queue_len=1"); // sanity
            queueP[gct].next_to_leave = buf_i;
        }

        // OOB: copy data to the buffer we added to the GetQBit queueP - dispatcher will consume it next, so queue
        // normally doesn't grow more than a handful of buffers
        if (gct == GCT_OOB) {
            buf_alloc (evb, &queueP[gct].gc_txts[buf_i], 0, segconf.vb_size, char, 0, NULL); // allocate to full size of VB so it doesn't need to be realloced
            buf_copy (evb, &queueP[gct].gc_txts[buf_i], &componentsP[comp_i].txt_data, char, 0, 0, 0);
        }

        // DEPN: data is only consumed after all MAIN and OOB data is consumed, so queue might grow very large.
        // We compress it so we can have more in memory and less offloaded to disk
        else {
            // copy data to another buffer, as this one will be freed and recycled
            buf_copy (evb, &depn.thread_data, &componentsP[comp_i].txt_data, uint8_t, 0, 0, NULL);

            int err;
            if (!is_final_flush && global_max_threads > 1) {
                // compress in separate thread, as to not delay releasing this VB. Only one thread is run in parallel as we join
                // the previous before creating a new one.
                ASSERT (!(err = pthread_create (&depn.thread, NULL, gencomp_compress_depn, &queueP[gct].gc_txts[buf_i])), 
                        "failed to create gencomp_compress_depn thread: %s", strerror(err));
                
                depn.has_thread = true;
            }
            else 
                gencomp_compress_depn (&queueP[gct].gc_txts[buf_i]);
        }
    }

    componentsP[comp_i].txt_data.len = componentsP[comp_i].txt_data.param = 0; // reset

    if (flag.debug_gencomp) debug_gencomp(comp_i==1 ? "flushcomp1" : "flushcomp2", false);
}

// called from compute_vb, with VBs in arbitrary order, and serializes them to their order. Notes:
// (1) We need to do it in the compute thread (rather than the main thread) so that zip_prepare_one_vb_for_dispatching, 
//     running in the main thread, can busy-wait for all MAIN compute threads to complete before flushing the final txt_data.
// (2) We do it in VB order, as recon_plan creation expects the lines to be in the same order as in the MAIN component. 
void gencomp_absorb_vb_gencomp_lines (VBlockP vb) 
{
    // skip VBs that are not MAIN, or have no gencomp lines
    if (vb->comp_i != COMP_MAIN || !vb->gencomp_lines.len) {        
        if (vb->comp_i == COMP_MAIN) num_MAIN_vbs_absorbedP++;
        
        serializer_skip (gc_vb_serializer, vb->vblock_i);
        return;
    } 

    serializer_lock (gc_vb_serializer, vb->vblock_i);

    mutex_lock (gc_protected);

    START_TIMER; // not including mutex wait times

    if (flag.debug_gencomp) iprintf ("absorb       %s/%u\n", comp_name (vb->comp_i), vb->vblock_i);

    ASSERT0 (!finished_absorbingP, "Absorbing is done - not expecting to be called");

    // limit vb_size of depn to 64MB, to reduce the gencomp write bottleneck, esp in --best in files with lots of depn
    // also: a bug in rans_compress_to_4x16 (called from compress_depn_buf) erroring in buffers near 2GB.
    uint32_t comp_size[3] = { [1] = segconf.vb_size,
                              [2] = componentsP[2].type == GCT_DEPN ? MIN_(segconf.vb_size, 1<<26) : segconf.vb_size };

    for (int i=1; i<=2; i++)
        buf_alloc (evb, &componentsP[i].txt_data, 0, comp_size[i], char, 1, "gencomp_txt_data");    

    // iterate on all lines the segmenter decided to send to gencomp (lines are of mixed componentsP)
    ARRAY (GencompLineIEntry, gc_lines, vb->gencomp_lines)
    for (uint32_t line_i=0; line_i < gc_lines_len; line_i++) { // cannot use for_buf here

        GencompLineIEntry *gcl = &gc_lines[line_i];

        // flush previous lines if there is no room for new line
        if (componentsP[gcl->comp_i].txt_data.len + gcl->line_len > comp_size[gcl->comp_i]) 
            gencomp_flush (gcl->comp_i, false);

        // note: it is possible that the some lines will fit into the queue, while subsequent
        // lines will be slated for re-reading 
        if (depn_method == DEPN_REREAD && 
            componentsP[gcl->comp_i].type == GCT_DEPN && 
            queueP[GCT_DEPN].next_unused == END_OF_LIST) {
            
            RereadLine reread_line = { .offset = gcl->offset, .line_len = gcl->line_len };
            buf_append_one (reread_depn_lines, reread_line);
            
            componentsP[GCT_DEPN].num_lines++; // for non-reread lines, incremented in gencomp_flush
        }

        // to to componentsP.txt_data to be flushed to the queue later
        else {
            buf_add (&componentsP[gcl->comp_i].txt_data, Btxt (gcl->line_index), gcl->line_len);
            componentsP[gcl->comp_i].txt_data.count++; 
        }
    }

    num_MAIN_vbs_absorbedP++;

    mutex_unlock (gc_protected);

    // we declare the file has having gencomp only if we actually have gencomp data 
    if (gc_lines_len)
        z_file->z_flags.has_gencomp = true; 

    COPY_TIMER (gencomp_absorb_vb_gencomp_lines);

    serializer_unlock (gc_vb_serializer);
}

//------------------------------------------------------------------------------------
// "Dispatcher functions" - Main thread functions called from within the dispatcher 
// loop - running in parallel with "Absorbing functions"
//------------------------------------------------------------------------------------

// main thread - called with gc_protected locked. remove the "next to leave" from the queue and place it
// on the unused stack. Caller should reset the buffer's len after copying it elsewhere. called with gc_protect locked.
static int gencomp_buf_leaves_queue (GencompType gct)
{
    int buf_i = queueP[gct].next_to_leave;

    // remove from buffer from queue
    queueP[gct].next_to_leave = GetQBit(buf_i, prev);
    if (queueP[gct].next_to_leave != END_OF_LIST) 
        SetQBit(queueP[gct].next_to_leave, next) = END_OF_LIST;
    else
        queueP[gct].last_added = END_OF_LIST; // buf_i was the head AND the tail of the list

    // add buffer to unused stack
    SetQBit(buf_i, next) = queueP[gct].next_unused;
    SetQBit(buf_i, prev) = END_OF_LIST; // always EOL for stack

    queueP[gct].next_unused = buf_i;

    return buf_i;
}

// main thread - called with gc_protected locked
static void gencomp_get_txt_data_from_queue (VBlockP vb, GencompType gct)
{
    int buf_i = gencomp_buf_leaves_queue (gct);
    BufferP buf = &queueP[gct].gc_txts[buf_i];

    if (gct == GCT_OOB)
        buf_copy (vb, &vb->txt_data, buf, char, 0, 0, "txt_data");
    else {
        // compressed buffer: first 4 bytes are uncomp_len, then the compressed data
        vb->txt_data.len = *B1ST32 (*buf);
        buf_alloc (vb, &vb->txt_data, 0, segconf.vb_size, char, 0, "txt_data"); 
    
        codec_rans_uncompress (evb, CODEC_RANS8, 0, buf->data + 4, buf->len - 4, 
                               &vb->txt_data, vb->txt_data.len, 0, "txt_data");
    }

    // initialize VB
    vb->comp_i    = GetQBit(buf_i, comp_i);
    vb->lines.len = GetQBit(buf_i, num_lines);    
    buf->len = 0;  // reset

    num_vbs_dispatched[gct]++;
    queueP[gct].queue_len--;

    if (flag.debug_gencomp) 
        debug_gencomp (vb->comp_i==1 ? "disp_comp1" : "disp_comp2", false);

    mutex_unlock (gc_protected);
} 

// main thread
static void gencomp_get_txt_data_from_disk (VBlockP vb)
{
    // case: we have data on disk - use it first
    if (!depn.offload_info.next)  // first read
        ASSERT (!fseek (depn.fp, 0, SEEK_SET), "fseek(%s,0) failed: %s", depn.name, strerror(errno));

    TxtDataInfoType *info = B(TxtDataInfoType, depn.offload_info, depn.offload_info.next++); 
    vb->comp_i         = info->comp_i;
    vb->lines.len      = info->num_lines;

    // copy data to VB
    buf_alloc (vb, &vb->txt_data, 0, segconf.vb_size, char, 0, "txt_data"); // allocate to full size of VB so it doesn't need to be realloced

    // note: thread_data_comp is already allocated, from when used for compression
    ASSERT (fread (depn.thread_data_comp.data, info->comp_len, 1, depn.fp) == 1, 
            "Failed to read buffer #%u length %u from %s", (int)depn.offload_info.next-1, (int)vb->txt_data.len, depn.name);
    depn.thread_data_comp.len = info->comp_len;

    vb->txt_data.len = *B1ST32 (depn.thread_data_comp); // first 4 bytes = uncomp_len
    ASSERT (vb->txt_data.len <= segconf.vb_size, "Invalid vb->txt_data.len=%"PRIu64, vb->txt_data.len); // sanity

    codec_rans_uncompress (evb, CODEC_RANS8, 0, depn.thread_data_comp.data+4, depn.thread_data_comp.len-4, 
                           &vb->txt_data, vb->txt_data.len, 0, "txt_data");

    if (flag.debug_gencomp) 
        iprintf ("Read from disk: buf=%"PRIu64" vb=%s num_lines=%u uncomp_len=%u comp_len=%u uncomp_alder32=%u comp_adler32=%u\n",
                  depn.offload_info.next-1, VB_NAME, info->num_lines, vb->txt_data.len32, info->comp_len, 
                  adler32 (1, STRb(vb->txt_data)), adler32 (1, STRb(depn.thread_data_comp))); 

    depn.thread_data_comp.len = depn.thread_data.len = 0;

    if (flag.debug_gencomp) debug_gencomp ("disp_comp2DSK", true);
}

// main thread - creates vb->reread_prescription with offsets of lines to be reread. the actual
// re-reading is done by the Seg compute thread for this vb
static void gencomp_prescribe_reread (VBlockP vb)
{
    ARRAY (RereadLine, lines, reread_depn_lines);

    buf_alloc (vb, &vb->reread_prescription, segconf.vb_size / lines[0].line_len, 0, RereadLine, 2, "reread_prescription"); // initial allocation, might grow

    vb->comp_i = SAM_COMP_DEPN;

    while (reread_depn_lines.next < lines_len) {
        if (vb->txt_data.len32 + lines[reread_depn_lines.next].line_len > segconf.vb_size) break; // VB is full

        buf_append_one (vb->reread_prescription, lines[reread_depn_lines.next]);
        vb->txt_data.len32 += lines[reread_depn_lines.next].line_len;
        vb->lines.len++;
        reread_depn_lines.next++;
    }

    if (flag.debug_gencomp) debug_gencomp ("disp_reread", true);
}

// main thread: populate vb->txt_data with the next buffer on the out-of-band or DEPN queue
bool gencomp_get_txt_data (VBlockP vb)
{
    #define DEBUG_GENCOMP(msg) \
        if (flag.debug_gencomp) \
            iprintf ("Returning %s vb=%s/%u with txt_data: len=%"PRIu64" adler32=%u\n", \
                    (msg), comp_name (vb->comp_i), vb->vblock_i, vb->txt_data.len, adler32 (1, STRb(vb->txt_data))); 

    #define RETURN(msg, call) \
        ( { call;\
            DEBUG_GENCOMP (msg);\
            zip_init_vb (vb);\
            return true; })

    if (!gc_protected.initialized) return false; // definitely no gencomp at all in this file (but if initialized, it doesn't always mean there is gencomp)

    mutex_lock (gc_protected);

    // case: we have no GetQBit data available, but MAIN data has all been dispatched (so no point returning
    // because there is no more file data to read), and some MAIN VBs have not sent data to absorb yet - 
    // 
    // case: we have out-of-band data: send it first
    if (queueP[GCT_OOB].queue_len) 
        RETURN ("OOB_Q", gencomp_get_txt_data_from_queue (vb, GCT_OOB)); // also unlocks mutex

    // case: finished ingesting PRIM, and finish all (if any) disk-offloaded data. Now we can do the in-memory GCT_DEPN queue   
    if (sam_finished_ingesting_prim && queueP[GCT_DEPN].queue_len)
        RETURN ("DEPN_Q", gencomp_get_txt_data_from_queue (vb, GCT_DEPN));

    // we might have data on disk - but we will not be accessing the queue or protected data anymore
    // and we don't want to hold the mutex while reading from disk
    mutex_unlock (gc_protected);

    // case: finished ingesting PRIM and no more out-of-band data, and all MAIN data has been flushed (which also means txt_file reached EOF,
    // see zip_prepare_one_vb_for_dispatching) - so no more MAIN or GetQBit data will be available. time for DEPN data.
    if (sam_finished_ingesting_prim && depn.offload_info.next < depn.offload_info.len) 
        RETURN ("DISK", gencomp_get_txt_data_from_disk (vb));

    // case: finished ingesting PRIM and no more out-of-band data, and all MAIN data has been flushed (which also means txt_file reached EOF,
    // see zip_prepare_one_vb_for_dispatching) - so no more MAIN or GetQBit data will be available. time for DEPN data.
    if (sam_finished_ingesting_prim && reread_depn_lines.next < reread_depn_lines.len) 
        RETURN (txt_file->codec == CODEC_BGZF ? "REREAD_BGZF" : "REREAD_PLAIN", gencomp_prescribe_reread (vb));

    // no more data exists at this point OR we have GCT_DEPN, but not finished ingesting PRIM yet
    // DEBUG_GENCOMP ("NO_DATA_AVAILABLE"); // commenting out because there are too many of
    return false; // no data available now

    #undef RETURN
}

// main thread: true if there is data on gencomp queues, or there might be in the future. called by
// zip_prepare_one_vb_for_dispatching after txt_file is exhausted.
bool gencomp_am_i_expecting_more_txt_data (void)
{
    if (!gc_protected.initialized) return false; // no gencomp at all in this file

    mutex_lock (gc_protected);

    if (!finished_absorbingP && num_vbs_dispatched[GCT_NONE] == num_MAIN_vbs_absorbedP) {
        // final flush. at this point, our thread is considered to be both "Dispatching" and "Absorbing"
        for (CompIType comp_i=1; comp_i <= MAX_GEN_COMP; comp_i++) 
            gencomp_flush (comp_i, z_sam_gencomp || comp_i==MAX_GEN_COMP); // final flush of gct, not of comp_i
        
        if (flag.debug_gencomp) iprintf ("Finished absorbing: num_MAIN_vbs_absorbedP=%u\n", num_MAIN_vbs_absorbedP);
        finished_absorbingP = true;
    }

    bool expecting = !finished_absorbingP || queueP[GCT_OOB].queue_len || queueP[GCT_DEPN].queue_len ||
                     reread_depn_lines.next < reread_depn_lines.len;

    if ((TXT_DT(DT_SAM) || TXT_DT(DT_BAM)) && finished_absorbingP && !queueP[GCT_OOB].queue_len && !num_vbs_dispatched[GCT_OOB]) {
        sam_finished_ingesting_prim = true;
        if (flag.debug_gencomp) iprint0 ("No PRIM VBs in this file\n");
    }

    mutex_unlock (gc_protected);
    return expecting;
} 

// main thread
void gencomp_a_main_vb_has_been_dispatched (void)
{
    if (!gc_protected.initialized) return; // no gencomp at all in this file

    num_vbs_dispatched[GCT_NONE]++;
}

// main thread - called be sam_zip_after_compute for PRIM VBs, in order of VBs. 
void gencomp_sam_prim_vb_has_been_ingested (VBlockP vb)
{
    num_SAM_PRIM_vbs_ingested++; // only accessed by this function    

    // note: we enforce the order between finished_absorbingP and queue_len in store, so we can be relaxed in load
    mutex_lock (gc_protected);
    bool my_finished_absorbing = finished_absorbingP; 
    uint16_t prim_queue_len = queueP[GCT_OOB].queue_len; // in SAM, GetQBit is PRIM
    mutex_unlock (gc_protected);

    if (flag.debug_gencomp) iprintf ("Ingested SA Groups of vb=%s\n", VB_NAME);

    // thread safety: 1. finished_absorbingP and  prim_queue_len these two are updated by the gencomp_absorb_vb_gencomp_lines 
    // which guarantees that if we have data, at least one of these two conditions will be true. 
    // 2. gencomp_get_txt_data ensures that if there is a VB being dispatched, it is accounted for in
    // at least one of: num_vbs_dispatched[GCT_OOB], queueP[GCT_OOB].queue_len
    if ((VB_DT(DT_SAM) || VB_DT(DT_BAM)) && my_finished_absorbing && !prim_queue_len && num_vbs_dispatched[GCT_OOB] == num_SAM_PRIM_vbs_ingested) {
        sam_sa_prim_finalize_ingest ();
        sam_finished_ingesting_prim = true;
        if (flag.debug_gencomp) iprintf ("Finished ingesting SA Groups: num_SAM_PRIM_vbs_ingested=%u\n", num_SAM_PRIM_vbs_ingested);
    }
}

// ZIP: compute thread of a DEPN VB: actually re-reading data into txt_data according to vb->reread_prescription
void gencomp_reread_lines_as_prescribed (VBlockP vb)
{
    buf_alloc_exact (vb, vb->txt_data, vb->txt_data.len, char, "txt_data");
    vb->txt_data.len = 0;

    // open a file handle private to this VB
    FILE *file = fopen (txt_file->name, "rb");
    ASSERT (file, "%s: Failed to open %s for rereading depn lines: %s", VB_NAME, txt_file->name, strerror(errno));

    if (txt_file->codec == CODEC_BGZF) 
        bgzf_reread_uncompress_vb_as_prescribed (vb, file);

    else { // CODEC_NONE
        for_buf (RereadLine, line, vb->reread_prescription) {
            ASSERT (!fseeko64 (file, line->offset.offset, SEEK_SET),
                    "%s: fseeko64 on %s failed while rereading depn lines: %s", VB_NAME, txt_file->name, strerror(errno));
            
            ASSERT (fread (BAFTtxt, line->line_len, 1, file) == 1,
                    "%s: fread of %u bytes on %s failed while rereading depn lines: %s", VB_NAME, line->line_len, txt_file->name, strerror(errno));

            vb->txt_data.len32 += line->line_len;
        }
    }

    fclose (file);
}