import os, sys

iteration = "Run4_2/"
datasets = ["RVSMPt10PU","RVSMPt10noPU","RVSMPt100PU","RVSMPt100noPU","RVSMPt1000PU","RVSMPt1000noPU","RVSMFlatPU","RVSMFlatnoPU"]
datasets = datasets + ["RVDMPt2PU","RVDMPt2noPU","RVDMPt10PU","RVDMPt10noPU","RVDMPt30PU","RVDMPt30noPU"]
# datasetsjoined = " ".join(datasets)

eosdir = "/eos/user/s/siluo/Muon/" + iteration
curdir = os.path.abspath(os.path.curdir)
subdir = os.path.join(curdir,"Submits")
dsdir = "/afs/cern.ch/work/s/siluo/Muon/filenames/"

# Make directories
if True:
  for ds in datasets:
    mklogdir = os.path.join(subdir,"logs"+ds)
    if not os.path.exists(mklogdir): os.makedirs(mklogdir)
    mkeosdir = eosdir+ds
    if not os.path.exists(mkeosdir): os.makedirs(mkeosdir)
    mkeosinprogressdir = os.path.join(mkeosdir,"InProgress")
    if not os.path.exists(mkeosinprogressdir): os.makedirs(mkeosinprogressdir)

# Make Submission scripts
if True:
  with open("SubTemplate.sub") as fin:
    lines = fin.readlines()
  for ds in datasets:
    fn = "Submits/"+ds+".sub"
    f = open(fn,"w")
    nf = len(open(dsdir+iteration+ds+".txt").readlines())
    writeline = lines
    writeline[1] = "arguments    = $(Proxy_path) $(ProcID) "+iteration+ds+" 0\n"
    logline = "logs"+ds+"/log_$(ClusterID)_$(ProcID)"
    writeline[4] = "output       = " + logline + ".out\n"
    writeline[5] = "error        = " + logline + ".err\n"
    writeline[6] = "log          = " + logline + ".log\n"
    writeline[12] = "queue "+str(nf)
    f.writelines(writeline)

if True:
  with open("MultiSub.sh") as fin:
    lines = fin.readlines()
  fn = "Submits/MultiSub.sh"
  f = open(fn,"w")
  lines[18] = "cd "+curdir+"/\n"
  f.writelines(lines)
  mode = "755"
  os.chmod(fn,int(mode,base=8));

if True:
  lines = []
  lines.append("#!/bin/sh\n")
  lines.append("for i in $(ls *.sub)\n")
  lines.append("do\n")
  lines.append("  echo Submitting $i\n")
  lines.append("  condor_submit $i\n")
  lines.append("done")
  fn = "Submits/Submission.sh"
  f = open(fn,"w")
  f.writelines(lines)
  mode = "755"
  os.chmod(fn,int(mode,base=8));
