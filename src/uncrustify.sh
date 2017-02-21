#! /bin/sh

files=$(find . -name '*.[ch]')

for item in $files ; do

    uncrustify -l C -c uncrustify.cfg $item
    mv $item.uncrustify $item

done
