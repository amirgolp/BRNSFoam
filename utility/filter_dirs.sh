#!/usr/bin/env bash
set -euo pipefail

########################################
# CONFIG
########################################

# Set to "true" for a dry run (only prints what would be deleted).
# Set to "false" to actually delete.
DRY_RUN=false

# Number of parallel delete workers
NUM_PROCS=8

# Directories to always keep
KEEP_ALWAYS=("0" "0_org" "constants" "system")

########################################
# HELPERS
########################################

is_always_keep() {
    local name="$1"
    for k in "${KEEP_ALWAYS[@]}"; do
        if [[ "$name" == "$k" ]]; then
            return 0
        fi
    done
    return 1
}

########################################
# MAIN
########################################

to_delete=()

for dir in */; do
    [ -d "$dir" ] || continue
    name="${dir%/}"

    # Always keep explicitly listed dirs
    if is_always_keep "$name"; then
        echo "KEEP (always): $name"
        continue
    fi

    action="keep"

    # Work on a copy without leading minus
    n="$name"
    [[ "$n" == -* ]] && n="${n#-}"

    # Scientific notation -> delete (implicit many decimals)
    if [[ "$n" == *[eE]* ]]; then
        action="delete"

    # Only digits and optional single dot
    elif [[ "$n" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
        if [[ "$n" == *.* ]]; then
            frac="${n#*.}"              # part after "."
            if ((${#frac} <= 3)); then
                # 0–4 decimal digits -> keep
                action="keep"
            else
                # Digits beyond the 4th decimal
                beyond="${frac:3}"
                if [[ "$beyond" =~ ^0+$ ]]; then
                    # only zeros beyond 4th -> effectively 4 decimals
                    action="keep"
                else
                    # non-zero in 5th/6th/... decimal place -> delete
                    action="delete"
                fi
            fi
        else
            # Integer (no decimal part) -> keep
            action="keep"
        fi
    else
        # Non-numeric names -> keep
        action="keep"
    fi

    if [[ "$action" == "keep" ]]; then
        echo "KEEP: $name"
    else
        echo "MARK FOR DELETE: $name"
        to_delete+=("$dir")
    fi
done

if ((${#to_delete[@]} == 0)); then
    echo "Nothing to delete."
    exit 0
fi

if $DRY_RUN; then
    echo
    echo "Dry run enabled (DRY_RUN=true). The following would be deleted:"
    printf '  %s\n' "${to_delete[@]}"
    exit 0
fi

echo
echo "Deleting ${#to_delete[@]} directories in parallel with $NUM_PROCS workers..."

# Parallel deletion
printf '%s\0' "${to_delete[@]}" | xargs -0 -P "$NUM_PROCS" -n 10 rm -rf --

echo "Done."
