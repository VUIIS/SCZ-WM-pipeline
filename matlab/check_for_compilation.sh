#/bin/bash
#
# Check if any of the matlab source code is newer than the compiled executable.

flag=0
files=`find src -name \* -type f`
for f in ${files}; do
    if [ "${f}" -nt bin/run_spm12.sh ]; then
        echo Source code newer than binary: "${f}"
        flag=1
    fi
done

if [ ${flag} == '1' ]; then
    exit 1
fi