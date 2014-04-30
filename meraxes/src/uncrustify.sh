#! /bin/sh

files=$(find . -name '*.[ch]')

for item in $files ; do

    uncrustify -l C -f $item -c uncrustify.cfg > $item.new
    mv $item.new $item

done
