#!/bin/bash
#
# Update the default namelist, streams, and stream_list files for a given
# MPAS core (namelist.CORE, streams.CORE, stream_list.CORE.*) in TRGT_DIR
# with the versions currently in SRCE_DIR.
#
# Any existing file in TRGT_DIR is first renamed to a time-stamped backup
# (e.g. namelist.CORE.YYYYMMDD_HHMMSS) before being replaced by the copy
# from SRCE_DIR.  If the new copy turns out to be identical to the backed-
# up old version, the backup is moved back into place instead, so that a
# file whose content hasn't actually changed keeps its original timestamp/
# inode rather than being needlessly replaced by an identical copy.
#
# Usage: update_default_inputs.sh SRCE_DIR TRGT_DIR CORE
#
set -eu
#
# Set the source and target directories.  Files (associated with the
# specified core) in the target directory will first be moved into time-
# stamped backup versions.  Then new files from the source directory will
# be copied from the source into the target directory.  Finally, this
# copying will be undone for only those files that did not need updating.
#
srce_dir="$1"
trgt_dir="$2"
#
# Get the MPAS core to be considered (e.g. init_atmosphere, atmosphere,
# etc).
#
core="$3"
#
# Set the extglob shell option, which enables extended pattern matching
# (e.g. +([!.])).  This is needed because a plain "*" glob can't express
# "one or more non-dot characters", so it would also match already time-
# stamped backup files (e.g. stream_list.CORE.SUFFIX.YYYYMMDD_HHMMSS) in
# addition to the plain files we actually want to match.
#
shopt -s extglob
#
# Get the time stamp to use when making back-ups of outdated files.
#

ts=$(date +%Y%m%d_%H%M%S)
#
# Loop through the namelist and stream files in the target directory
# (excluding any files that have time stamps in their names).   Rename
# each file by adding a time stamp at the end of its name.
#
cd "${trgt_dir}"
for file in namelist.${core} streams.${core} stream_list.${core}.+([!.]); do
    # Check for file existence because bash's behavior, when a list item
    # is a literal name or an unmatched glob/extglob pattern, is to leave
    # it in the loop as-is rather than dropping it -- so without this
    # check, mv below would be attempted on a nonexistent file (e.g. on a
    # from-scratch build where trgt_dir doesn't have these files yet),
    # which would abort the script since "set -e" is in effect.
    if [ -f "${file}" ]; then
        file_old="${file}.${ts}"
        mv "${file}" "${file_old}"
    fi
done
#
# Loop through the namelist and stream files in the source directory
# (excluding any files that have time stamps in their names) and:
#
# 1) Copy each file into the target directory.  This copy is at file_new_fp,
#    and we refer to it as the "new" file.
# 2) Check if a corresponding time-stamped file (i.e. an old version of
#    the file) was created in this directory in the previous loop.  If
#    so, compare the old file to the new file.  If they are identical,
#    overwrite the new file with the old one (and delete the old) since
#    there was no reason to update the old file.  If they are not identical,
#    then the file needed updating, so keep the new file and the time-
#    stamped old file as they are.
#
cd "${srce_dir}"
for file in namelist.${core} streams.${core} stream_list.${core}.+([!.]); do
    # Check for file existence for the same reason as in the loop above:
    # a literal name or unmatched glob/extglob pattern stays in the list
    # as-is, so without this check, cp below would be attempted on a
    # nonexistent file and abort the script (due to "set -e").
    if [ -f "${file}" ]; then
        file_new_fp="${trgt_dir}/${file}"
        cp "${file}" "${file_new_fp}"
        file_old_fp="${trgt_dir}/${file}.${ts}"
        if [ -f "${file_old_fp}" ]; then
            if diff "${file_new_fp}" "${file_old_fp}" > /dev/null; then
                mv "${file_old_fp}" "${file_new_fp}"
            fi
        fi
    fi
done

