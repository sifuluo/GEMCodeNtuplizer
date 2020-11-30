#!/bin/sh
export X509_USER_PROXY=$1
voms-proxy-info -all
voms-proxy-info -all -file $1
j=${4:-0}
i=$(($2+5000*$j))

existname=/eos/user/s/siluo/Muon/${3}/out_${i}.root
if [ -e $existname ]
then
  echo 'File already processed'
  exit 0
fi

# dir=$(pwd)
cd /afs/cern.ch/work/s/siluo/CMSSW/CMSSW_11_2_0_pre7/src
export SCRAM_ARCH=slc7_amd64_gcc700
eval `scramv1 runtime -sh`
cd /afs/cern.ch/work/s/siluo/Muon/Nov30Batch/
cmsRun NtuplizeSIM.py ifile=$i dataset=$3

filename=/eos/user/s/siluo/Muon/${3}/InProgress/out_${i}.root
filesize=$(stat -c %s $filename)
if [ $filesize -gt 1000 ]
then
  mv $filename $existname
  echo 'Job is done...'
else
  rm $filename
  echo "Job failed...Output filesize $filesize. Now retrying......"
  exit 1
fi
