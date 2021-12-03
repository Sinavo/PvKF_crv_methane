#!/bin/bash
####################################################################
############## LOOP TO RUN FOR SEVERAL DAYS/WEEKS ##################
####################################################################
#export IOAPI_DIR=/cm/shared/apps/ioapi/pgi/2015/64/3.1/Linux2-x86-64_pgi2015_netcdf3.6.3/

#for ((i=3;i<6;i+=1));do
echo 'i=' $i
###################################################
######## ALL INPUTS ARE SET HERE ##################
###################################################
export INITIALPDATE=20100401  #First day
export INITIALSDATE=20100402  #First day 
export INITIALEDATE=20100403  #Second day 
#export RUNDAYS=2
#export RUN=daily #daily, weekly, or 35 change as approperiate

export SSDATE=`datshift $INITIALSDATE $i | tail -5 | head -1 | cut -d " " -f1`
export SHH=00                                           #Start Hour

export EDATE=`datshift $INITIALEDATE $i | tail -5 | head -1 | cut -d " " -f1`
export EHH=00                                   #End Hour

export PDATE=`datshift $INITIALPDATE $i | tail -5 | head -1 | cut -d " " -f1`
export PHH=00
###################################################
### CUTTING DATE TO YEAR, MONTH, DAY, and HOUR ####
###################################################
export SMM=`echo $SSDATE | cut -c5-6`
export SDD=`echo $SSDATE | cut -c7-8`
export SYY=`echo $SSDATE | cut -c1-4`

echo $SSDATE $SHH $SMM $SDD $SYY
export EMM=`echo $EDATE | cut -c5-6`
export EDD=`echo $EDATE | cut -c7-8`
export EYY=`echo $EDATE | cut -c1-4`

export PMM=`echo $PDATE | cut -c5-6`
export PDD=`echo $PDATE | cut -c7-8`
export PYY=`echo $PDATE | cut -c1-4`

#export SHELL_RUN=/home/sinavo/pkg/cmaq/5.3/CCTM/scripts_var/var_1
export SHELL_RUN=/home/sinavo/pkg/cmaq/5.3/CCTM/da_9/var
###################################################
######## RUN HOUR ##################
###################################################
#      for j in {0..23..1}; do
          echo 'j=' $j
          if (( 10#$j < 10 )); then
          hour=0${j}0000
          else
          hour=${j}0000
          fi
echo "hour=" $hour
          b2=$(( 10#$j-1 ))
          if (( $b2 < 10 )); then
              hourb2=0${b2}0000
          else
              hourb2=${b2}0000
          fi
echo "hourb2=" $hourb2
#          echo ADDL = $hour
#          echo JTIME=ADDL0= $hourb2
          export COUNTER=$z
          export ADDL=$hour
          export ADDL0=$hourb2
         # export JTIME=$hourb2
          export JTIME=$hour
          export ADDL=`echo "$ADDL" | tail -5 | head -1 | cut -d " " -f1`
          export ADDL0=`echo "$ADDL0" | tail -5 | head -1 | cut -d " " -f1`
          export JTIME=`echo "$JTIME" | tail -4 | head -1 | cut -d " " -f1`
#          export COUNTER=`echo "$COUNTER" | tail -5 | head -1 | cut -d " " -f1`
          echo $ADDL $ADDL0 $APPL $APPL0
          sed "s:HRNOW:$ADDL:g" $SHELL_RUN/run_hr_hemi.csh | sed "s:HRPAST:$ADDL0:g" | sed "s:SYY:$SYY:g" | sed "s:SMM:$SMM:g" | sed "s:SDD:$SDD:g" | sed "s:PYY:$PYY:g" | sed "s:PMM:$PMM:g" | sed "s:PDD:$PDD:g"  > $SHELL_RUN/run.hr

         chmod 755 $SHELL_RUN/run.hr
         $SHELL_RUN/run.hr
#       done
###################################################
################### RUN DAYS ######################
###################################################
#export V1=`echo $SDATE | cut -c6-8`
#export V2=`echo $EDATE | cut -c6-8`
#sed "s:VAPPL:$V1$V2:g" $SHELL_RUN/run_cctm_hemi.csh | sed "s:SYY:$SYY:g" | sed "s:SMM:$SMM:g" | sed "s:SDD:$SDD:g" | sed "s:EYY:$EYY:g" | sed "s:EMM:$EMM:g" | sed "s:EDD:$EDD:g" | sed "s:SMM:$PYY:g" | sed "s:PMM:$PMM:g" | sed "s:PDD:$PDD:g"> run.cctm

#chmod 755 ./run.cctm
#./run.cctm

#done
