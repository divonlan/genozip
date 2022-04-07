#!/usr/bin/bash 

# ------------------------------------------------------------------
#   dict_id_gen.sh
#   Copyright (C) 2021-2022 Genozip Limited
#   Please see terms and conditions in the file LICENSE.txt

# dict_id_gen.h generation:
# Step 1: dict_id_gen.sh generates dict_id_gen.c, including all the GENDICT definitions from the data type include files (eg vcf.h)
# Step 2: dict_id_gen.sh compiles dict_id_gen.c and generate dict_id_gen[.exe]
# Step 3. dict_id_gen.sh generates dict_id_gen.h: it uses dict_id_gen[.exe]to generate the field constant, and then adds the fields enum and mapping

# The source information for this script are the GENDICT and GENDICT_PREFIX pragmas in the code

generate_dict_id_gen_c() 
{ 
cat > dict_id_gen.c << END
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <unistd.h> 

typedef enum { DTYPE_FIELD, DTYPE_1, DTYPE_2 } DictIdType;

#pragma pack(1) // structures that are part of the genozip format are packed.
#define DICT_ID_LEN 8
typedef union DictId { uint64_t num; uint8_t id[DICT_ID_LEN]; } DictId;
#pragma pack()

DictId dict_id_make (const char *str, unsigned str_len, DictIdType dict_id_type) 
{ 
    DictId dict_id = {0}; 

    if (!str_len) str_len = strlen (str);

    if (str_len <= DICT_ID_LEN) 
        memcpy (dict_id.id, str, str_len);
    
    else { 
        #define half1_len (DICT_ID_LEN/2)
        #define half2_len (DICT_ID_LEN - DICT_ID_LEN/2)

        memcpy (dict_id.id, str, half1_len); 
        memcpy (dict_id.id + half1_len, str+str_len-half2_len, half2_len);
    }

    switch (dict_id_type) {
        case DTYPE_FIELD  : dict_id.id[0] = dict_id.id[0] & 0x3f; break;
        case DTYPE_1      : dict_id.id[0] = dict_id.id[0] | 0xc0; break;
        case DTYPE_2      : break;
        default: fprintf (stderr, "Error in dict_id_gen.c/dict_id_make: invalid type %d", dict_id_type);
    }
    return dict_id;
}

int main (int argc, const char **argv) 
{
END

    max_fields=0 

    #pragma GENDICT REF_CONTIG=DTYPE_FIELD=CONTIG 
    tags=`egrep -w "^#pragma GENDICT" ${files[*]} | cut -d" " -f3 | tr -d '\15' ` # remove windows \r
    num_tags=`echo $tags | wc -w`

    for t in $tags; do

        IFS='=' # split separator (much, much faster than invoking cut)
        read -a split <<< "$t" 
        
        var=${split[0]}
        type=${split[1]}
        name=${split[2]}

        printf "    printf (\"#define _%s ((uint64_t)%%\"PRId64\")\\\\n\", dict_id_make (\"%s\", %s, %s).num);\n" $var $name ${#name} $type >> dict_id_gen.c
    done

    echo "}" >> dict_id_gen.c
} # generate_dict_id_gen_c

generate_dict_id_gen_h() 
{ 
    cat > dict_id_gen.h << END
// ------------------------------------------------------------------
// dict_id_gen.h
//   Copyright (C) 2019-2022 Genozip Limited
//   Please see terms and conditions in the file LICENSE.txt

// THIS FILE IS AUTO-GENERATED BY dict_id_gen.sh

#pragma once

#define TAG(tag) #tag, sizeof #tag - 1

END

    # add field constant definitions to dict_id_gen.h
    ./$dict_id_gen_exe >> dict_id_gen.h
    echo >> dict_id_gen.h

    max_fields=0 

    # add MAPPING
    IFS='=' # split separator
    for f in ${files[@]} ; do

        prefix=`egrep -w "^#pragma GENDICT_PREFIX" $f | head -1 | cut -d" " -f3 | tr -d '\15'`
        if [ ${#prefix} -eq 0 ]; then echo "dict_id_gen.sh: Error: no GENDICT_PREFIX in $f"; exit 1; fi

        # note: readarray doesn't work on MacOS :(
        readarray -t vars <<< `egrep -w "^#pragma GENDICT" $f | cut -d" " -f3 | cut -d"=" -f1`
        readarray -t tags <<< `egrep -w "^#pragma GENDICT" $f | cut -d" " -f3 | cut -d"=" -f3 | cut -d"/" -f1 | sed 's/[^\!-~]//g'`

        # calculate MAX_NUM_FIELDS_PER_DATA_TYPE
        if [ ${#vars[@]} -gt $max_fields ]; then max_fields=${#vars[@]}; fi

        # add enum: typedef enum { REF_CONTIG, ... , NUM_REF_FIELDS } REFFields;
        # add: typedef enum { REF_CONTIG, ... } Fields;
        printf "typedef enum { " >> dict_id_gen.h
        for v in ${vars[*]}; do
            echo -n "$v, " >> dict_id_gen.h
        done
        printf "NUM_%s_FIELDS } %sFields;\n\n" $prefix $prefix >> dict_id_gen.h
        
        # add MAPPING: [did_i]={ .num = _##did_i }
        echo "#define ${prefix}_PREDEFINED { \\"  >> dict_id_gen.h

        for ((i = 0 ; i < ${#vars[@]} ; i++)); do
            echo "    [${vars[$i]}] = { { _${vars[$i]} }, TAG(${tags[$i]}) }, \\" >> dict_id_gen.h
        done

        printf "} \n\n" >> dict_id_gen.h
    done

    # add: MAX_NUM_FIELDS_PER_DATA_TYPE
    printf "#define MAX_NUM_FIELDS_PER_DATA_TYPE %u\n\n" $max_fields >> dict_id_gen.h
}

is_windows=`uname|grep -i mingw`
if [ -n "$is_windows" ]; then
    dict_id_gen_exe=dict_id_gen.exe
else
    dict_id_gen_exe=dict_id_gen
fi

headers=(`egrep "^#include" data_types.h | cut -d\" -f2`)
files=(`egrep "^#pragma GENDICT" ${headers[*]} | cut -d: -f1 | uniq`)

generate_dict_id_gen_c

# We accept CC as $1, used in conda when called from Makefile
CC=$1
if [ ${#CC} -eq 0 ]; then CC=/mingw64/bin/gcc; fi
    
$CC dict_id_gen.c -o $dict_id_gen_exe

generate_dict_id_gen_h 

rm -f dict_id_gen.c
