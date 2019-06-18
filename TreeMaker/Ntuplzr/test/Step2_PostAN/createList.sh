#!/bin/bash

if [ "$#" -ne 4 ] ; then
  echo "Usage: $0 Sample V2/ext user" >&2
 # exit 1
fi

sample=$1
ext=$2
name=$3
iter=$4

dir=`ls /eos/cms/store/group/upgrade/RTB/${sample}/crab_${sample}/| awk '{print $1}'`

echo $sample
echo $dir
echo $ext
echo $name
ls /eos/cms/store/group/upgrade/RTB/${sample}/crab_${sample}/${dir}/0000/ | grep -v log | grep -v failed | awk -v dir=$dir -v sample=$sample -v ext=$ext -v name=$name -v iter=$iter '{print "/eos/cms/store/group/upgrade/RTB/" sample "/crab_" sample "/" dir"/0000/"$1}' >> input_${sample}.list 



