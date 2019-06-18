echo $1
root -l -b -q runAll.C\(\"input_$1.list\",\"$1\",$2\);
#root -l runAll.C\(\"input_$1.list\",\"$1\",$2\);
