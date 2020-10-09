#!/usr/bin/bash

WORK_JEDI=$1

cp fv3-bundle/fv3-jedi/test/CMakeLists.txt $WORK_JEDI/src/fv3-bundle/fv3-jedi/test/CMakeLists.txt
for FILE in convertstate_gfs_c96 dirac_gfs_bal_rf_c96 dirac_gfs_bal_mgbf_c96 
do 
  cp fv3-bundle/fv3-jedi/test/testinput/${FILE}.yaml $WORK_JEDI/src/fv3-bundle/fv3-jedi/test/testinput/${FILE}.yaml
  cp fv3-bundle/fv3-jedi/test/testoutput/${FILE}.ref $WORK_JEDI/src/fv3-bundle/fv3-jedi/test/testout/${FILE}.ref
done

cp fv3-bundle/saber/src/saber/CMakeLists.txt $WORK_JEDI/src/fv3-bundle/saber/src/saber/CMakeLists.txt
cp -R fv3-bundle/saber/src/saber/mgbf $WORK_JEDI/src/fv3-bundle/saber/src/saber/mgbf
for FILE in CMakeLists.txt ErrorCovarianceMGBF.h ParametersMGBF.h OoMgbf.h instantiateCovarFactory.h
do 
  cp fv3-bundle/saber/src/saber/oops/$FILE $WORK_JEDI/src/fv3-bundle/saber/src/saber/oops/$FILE
done

exit 0
