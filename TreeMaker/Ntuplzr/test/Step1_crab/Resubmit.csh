#!bin/tcsh
foreach i (`ls | grep crab_`)
  echo $i
  crab resubmit -d $i  
  #crab status $i  
end
