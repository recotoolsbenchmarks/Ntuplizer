#!bin/tcsh
foreach i (`ls | grep crab_`)
  echo $i
  #crab resumit $i  
  crab status $i  
end
