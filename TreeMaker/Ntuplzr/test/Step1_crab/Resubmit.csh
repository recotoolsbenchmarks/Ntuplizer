#!bin/tcsh
foreach i (`ls | grep crab_`)
  echo $i
  crab status $i  
  crab resubmit -d $i --maxmemory=2000 
end
