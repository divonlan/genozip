#!/usr/bin/env bash

# ------------------------------------------------------------------
#   generate_context_validate.sh
#   Copyright (C) 2026-2026 Genozip Limited. Patent Pending.
#   Please see terms and conditions in the file LICENSE.txt
#
#   WARNING: Genozip is proprietary, not open source software. Modifying
#   the source code is strictly prohibited and subject to penalties
#   specified in the license.
# ------------------------------------------------------------------

if [[ "$GENOZIP_HOME" == "" ]]; then # Windows note: definition in /home/divon/.bashrc overrides definition in Windows Settings->Environment Variables
    GENOZIP_HOME=~/genozip
fi

set -euo pipefail

HEADER="${GENOZIP_HOME}/src/context_struct.h"

[[ -f "$HEADER" ]] || exit 1

cat <<'EOF'
// ------------------------------------------------------------------
//   context_validate.c
//   Copyright (C) 2026-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying
//   the source code is strictly prohibited and subject to penalties
//   specified in the license.
// ------------------------------------------------------------------

#include <stddef.h>
#include "context_struct.h"

#define ASSERT_CONTEXT_OFFSET(field, expected_word) \
    ASSERT (offsetof(Context, field) == 8 * (expected_word), \
            "Context " #field ": expected word %zu (byte %zu), " \
            "actual word %zu (byte %zu)", \
            (size_t)(expected_word), \
            8 * (size_t)(expected_word), \
            offsetof(Context, field) / 8, \
            offsetof(Context, field))

#define ASSERT_BUFFER_ALIGNED_64(field) \
    ASSERT (offsetof(Context, field) % 64 == 0, \
            "Context Buffer " #field ": byte offset %zu (word %zu) " \
            "is not 64-byte aligned", \
            offsetof(Context, field), \
            offsetof(Context, field) / 8)

void context_validate (void)
{
EOF

awk '
function strip_comment(s,    p) {
    p = index(s, "//")
    if (p)
        s = substr(s, 1, p - 1)
    return s
}

function trim(s) {
    sub(/^[[:space:]]+/, "", s)
    sub(/[[:space:]]+$/, "", s)
    return s
}

# Is this a declaration that can consume a pending § marker?
#
# We deliberately accept any statement ending in ; here. This also
# recognizes declaration macros such as:
#
#     STR (last_snip);
#
# but does not recognize union/struct/enum opening braces.
function is_declaration(s) {
    s = trim(strip_comment(s))

    if (s == "")
        return 0

    if (s ~ /[{}]/)
        return 0

    return s ~ /;[[:space:]]*$/
}

# Extract the first declarator from a normal C declaration.
#
# Examples:
#
#   uint32_t a, b;       -> a
#   uint8_t x;           -> x
#   Codec a;             -> a
#   Buffer foo;          -> foo
#
# This intentionally only extracts the FIRST declarator.
function first_field(s,    t, semi, decl, n, a, i, token) {
    s = trim(strip_comment(s))

    sub(/;[[:space:]]*$/, "", s)

    # STR(...) is handled separately.
    if (s ~ /^[[:space:]]*STR[[:space:]]*\(/)
        return ""

    # Remove an initializer, if any.
    sub(/[[:space:]]*=[^,]*$/, "", s)

    # First declarator is before the first comma.
    semi = index(s, ",")
    if (semi)
        decl = substr(s, 1, semi - 1)
    else
        decl = s

    decl = trim(decl)

    # A bitfield is deliberately not passed to offsetof().
    if (decl ~ /:[[:space:]]*[0-9]+[[:space:]]*$/)
        return ""

    # Find the final identifier in the first declarator.
    n = split(decl, a, /[[:space:]]+/)

    for (i = n; i >= 1; i--) {
        token = a[i]

        # Remove pointer/array syntax around the identifier.
        gsub(/^[*]+/, "", token)
        sub(/\[[^]]*\]$/, "", token)

        if (token ~ /^[A-Za-z_][A-Za-z0-9_]*$/)
            return token
    }

    return ""
}

# Recognize bitfield declarations so that the marker is consumed without
# generating an invalid offsetof() expression.
function is_bitfield(s) {
    s = trim(strip_comment(s))
    sub(/;[[:space:]]*$/, "", s)
    return s ~ /:[[:space:]]*[0-9]+[[:space:]]*$/
}

BEGIN {
    in_context = 0
    depth = 0
    pending_word = -1
}

{
    raw = $0
    code = strip_comment(raw)

    # Find the beginning of typedef struct Context.
    if (!in_context) {
        if (code ~ /^[[:space:]]*typedef[[:space:]]+struct[[:space:]]+Context[[:space:]]*\{/) {
            in_context = 1
            depth = 1
        }
        next
    }

    # A § marker can occur in a comment, including the special typo §s.
    #
    # We use only the FIRST number after §.
    if (match(raw, /§[s]?[[:space:]]*[0-9]+/)) {
        marker = substr(raw, RSTART, RLENGTH)
        sub(/^.*§[s]?[[:space:]]*/, "", marker)
        pending_word = marker + 0
    }

    # Consume the pending marker on the first declaration that follows it.
    if (pending_word >= 0 && is_declaration(code)) {
        declaration = trim(strip_comment(code))

        # STR(last_snip); expands to a declaration and therefore does not
        # provide a field that can safely be named in offsetof(Context,...).
        if (declaration ~ /^STR[[:space:]]*\(/) {
            print "    // §" pending_word ": STR declaration macro -- offsetof check unavailable"
            pending_word = -1
        }

        # Bitfields cannot be operands of offsetof().
        else if (is_bitfield(declaration)) {
            print "    // §" pending_word ": bit-field -- offsetof check unavailable"
            pending_word = -1
        }

        else {
            field = first_field(declaration)

            if (field != "") {
                print "    ASSERT_CONTEXT_OFFSET(" field ", " pending_word ");"
                pending_word = -1
            }
        }
    }

    # Keep track of braces so we know when typedef struct Context ends.
    opens = gsub(/\{/, "{", code)
    closes = gsub(/\}/, "}", code)

    depth += opens - closes

    if (depth <= 0)
        exit
}

' "$HEADER"

cat <<'EOF'

    /* Every Buffer in Context must begin on a 64-byte boundary. */
EOF

awk '
function strip_comment(s,    p) {
    p = index(s, "//")
    if (p)
        s = substr(s, 1, p - 1)
    return s
}

BEGIN {
    in_context = 0
    depth = 0
}

{
    raw = $0
    code = strip_comment(raw)

    if (!in_context) {
        if (code ~ /^[[:space:]]*typedef[[:space:]]+struct[[:space:]]+Context[[:space:]]*\{/) {
            in_context = 1
            depth = 1
        }
        next
    }

    # Match Buffer declarations at any nesting level inside Context.
    if (code ~ /^[[:space:]]*Buffer[[:space:]]+[A-Za-z_][A-Za-z0-9_]*/) {
        line = code

        sub(/^[[:space:]]*Buffer[[:space:]]+/, "", line)
        sub(/[[:space:]]*;.*/, "", line)

        # Buffer declarations in Context are normally one declarator.
        # Keep only the first one if there is ever a comma-separated
        # declaration.
        comma = index(line, ",")
        if (comma)
            line = substr(line, 1, comma - 1)

        sub(/\[[^]]*\]$/, "", line)

        print "    ASSERT_BUFFER_ALIGNED_64(" line ");"
    }

    opens = gsub(/\{/, "{", code)
    closes = gsub(/\}/, "}", code)

    depth += opens - closes

    if (depth <= 0)
        exit
}

' "$HEADER"

cat <<'EOF'

    /* Verify sizes of project types actually used as Context fields. */
    ASSERT (sizeof (Buffer) == 64, "sizeof(Buffer)=%zu, expected 64", sizeof (Buffer));

    ASSERT (sizeof (DictId) == 8, "sizeof(DictId)=%zu, expected 8", sizeof (DictId));
    ASSERT (sizeof (Did) == 2, "sizeof(Did)=%zu, expected 2", sizeof (Did));
    ASSERT (sizeof (LocalType) == 1, "sizeof(LocalType)=%zu, expected 1", sizeof (LocalType));
    ASSERT (sizeof (struct FlagsCtx) == 1, "sizeof(struct FlagsCtx)=%zu, expected 1", sizeof (struct FlagsCtx));
    ASSERT (sizeof (struct FlagsDict) == 1, "sizeof(struct FlagsDict)=%zu, expected 1", sizeof (struct FlagsDict));
    ASSERT (sizeof (B250Size) == 1, "sizeof(B250Size)=%zu, expected 1", sizeof (B250Size));
    ASSERT (sizeof (Codec) == 1, "sizeof(Codec)=%zu, expected 1", sizeof (Codec));

    ASSERT (sizeof (ValueType) == 8, "sizeof(ValueType)=%zu, expected 8", sizeof (ValueType));
    ASSERT (sizeof (WordIndex) == 4, "sizeof(WordIndex)=%zu, expected 4", sizeof (WordIndex));
    ASSERT (sizeof (TxtWord) == 8, "sizeof(TxtWord)=%zu, expected 8", sizeof (TxtWord));
    ASSERT (sizeof (LineIType) == 4, "sizeof(LineIType)=%zu, expected 4", sizeof (LineIType));
    ASSERT (sizeof (IdType) == 1, "sizeof(IdType)=%zu, expected 1", sizeof (IdType));
    ASSERT (sizeof (SamFlags) == 2, "sizeof(SamFlags)=%zu, expected 2", sizeof (SamFlags));
    ASSERT (sizeof (PosType32) == 4, "sizeof(PosType32)=%zu, expected 4", sizeof (PosType32));
    ASSERT (sizeof (thool) == 1, "sizeof(thool)=%zu, expected 1", sizeof (thool));
    ASSERT (sizeof (Ploidy) == 1, "sizeof(Ploidy)=%zu, expected 1", sizeof (Ploidy));
    ASSERT (sizeof (ContextP) == 8, "sizeof(ContextP)=%zu, expected 8", sizeof (ContextP));
    ASSERT (sizeof (SnipIterator) == 8, "sizeof(SnipIterator)=%zu, expected 8", sizeof (SnipIterator));
    ASSERT (sizeof (BamAssTrimCigarTreatment) == 1, "sizeof(BamAssTrimCigarTreatment)=%zu, expected 1", sizeof (BamAssTrimCigarTreatment));
    ASSERT (sizeof (StoreType) == 1, "sizeof(StoreType)=%zu, expected 1", sizeof (StoreType));
    ASSERT (sizeof (LocalDepType) == 1, "sizeof(LocalDepType)=%zu, expected 1", sizeof (LocalDepType));
    ASSERT (sizeof (ContainerP) == 8, "sizeof(ContainerP)=%zu, expected 8", sizeof (ContainerP));
    ASSERT (sizeof (SectionType) == 1, "sizeof(SectionType)=%zu, expected 1", sizeof (SectionType));
    ASSERT (sizeof (SpecialResult) == 1, "sizeof(SpecialResult)=%zu, expected 1", sizeof (SpecialResult));

    /* Types defined inside Context. */
    ASSERT (sizeof (struct ctx_tp) == 8,
            "sizeof(struct ctx_tp)=%zu, expected 8", sizeof (struct ctx_tp));

    ASSERT (sizeof (Context) % 64 == 0,
            "sizeof(Context)=%zu (8B-words=%1.3f) is not a multiple of 64: %u %% 64 = %u",
            sizeof (Context), (double)sizeof (Context) / 8.0, sizeof (Context), sizeof (Context) % 64);            

    ASSERT (RECON_STATE_SIZE <= 61, "Expecting RECON_STATE_SIZE=%u <= 61 (see Ice struct)", RECON_STATE_SIZE);

    if (flag.show_memory)
        iprintf ("Context-related sizes: Context=%u (%u x 64B) Mutex=%u RECON_STATE_SIZE=%u\n", 
                 (unsigned)sizeof(Context), (unsigned)sizeof(Context)/64, (unsigned)sizeof(Mutex), (unsigned)RECON_STATE_SIZE);
}
EOF
