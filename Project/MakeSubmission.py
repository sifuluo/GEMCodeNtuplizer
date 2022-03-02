import os, sys

InPath = "Private/"
OutPath = "MuGun/Ntuple/"
cmssw_ver = "CMSSW_12_2_0_pre3"
maxfile   = 5000
#Run4
# datasets = ["RVSMPt10PU","RVSMPt10noPU","RVSMPt100PU","RVSMPt100noPU","RVSMPt1000PU","RVSMPt1000noPU","RVSMFlatPU","RVSMFlatnoPU"]
# datasets = datasets + ["RVDMPt2PU","RVDMPt2noPU","RVDMPt10PU","RVDMPt10noPU","RVDMPt30PU","RVDMPt30noPU"]

#Run3
# datasets = ["DY","H9","H15","H30","H100"]

#Private
# datasets = ["SMNoPU0","SMNoPU0a","SMNoPU1","SMNoPU1a","SMNoPU2","SMNoPU2a"]
datasets = ["SMNoPU0","SMNoPU1","SMNoPU2","DMNoPU0","DMNoPU1","DMNoPU2"]

if len(sys.argv) == 1:
  process = 0
else:
  process = int(sys.argv[1])

eosdir = "/eos/user/s/siluo/" + OutPath
curdir = os.path.abspath(os.path.curdir)
subdir = os.path.join(curdir,"Submits")
dsdir = "/afs/cern.ch/work/s/siluo/Muon/filenames/"

# Make directories
if process == 0 or process == 1:
  for ds in datasets:
    mklogdir = os.path.join(subdir,"logs"+ds)
    if not os.path.exists(mklogdir):
      os.makedirs(mklogdir)
      print('  Created: ' + ds + ' log directory: ' + mklogdir)
    else: print( '  Existed: ' + ds + ' log directory: ' + mklogdir)

    mkeosdir = eosdir+ds
    if not os.path.exists(mkeosdir):
      os.makedirs(mkeosdir)
      print('  Created: ' + ds + ' eos directory: ' + mkeosdir)
    else: print( '  Existed: ' + ds + ' eos directory: ' + mkeosdir)

    mkeosinprogressdir = os.path.join(mkeosdir,"InProgress")
    if not os.path.exists(mkeosinprogressdir):
      os.makedirs(mkeosinprogressdir)
      print('  Created: ' + ds + ' WIP directory: ' + mkeosinprogressdir)
    else: print( '  Existed: ' + ds + ' WIP directory: ' + mkeosinprogressdir)

    print('Dataset: ' + ds + ' finished creating directories')
  print('Iteration: ' + OutPath + ' finished creating directories')

# Make Submission scripts
if process == 0 or process == 2:
  for ds in datasets:
    nf = len(open(dsdir+InPath+ds+".txt").readlines())
    lines = []
    lines.append("Proxy_path   = /afs/cern.ch/user/s/siluo/x509up\n")
    lines.append("arguments    = $(Proxy_path) $(ProcID) "+InPath+ds+" 0 "+OutPath+ds+"\n")
    lines.append("executable   = MultiSub.sh\n")
    lines.append("max_retries  = 10\n")
    lines.append("+JobBatchName= " + InPath+ds +"\n")
    lines.append("output       = logs"+ds+"/log_$(ClusterID)_$(ProcID).out\n")
    lines.append("error        = logs"+ds+"/log_$(ClusterID)_$(ProcID).err\n")
    lines.append("log          = logs"+ds+"/log_$(ClusterID)_$(ProcID).log\n")
    lines.append("universe     = vanilla\n")
    lines.append('Requirements = (OpSysAndVer =?= "CentOS7")\n')
    lines.append('+JobFlavour  = "longlunch"\n')
    lines.append("RequestCpus  = 2\n")
    lines.append("periodic_release =  (NumJobStarts < 10) && ((CurrentTime - EnteredCurrentStatus) > (5*60))\n")
    lines.append("queue "+str(nf)+"\n")

    fn = "Submits/"+ds+".sub"
    if nf > maxfile:
      fn = "Submits/"+ds+"0.sub"
      lines[13] = "queue "+str(maxfile)+"\n"
    f = open(fn,"w")
    f.writelines(lines)
    print("Dataset: " + ds + ' submission script created: '+ fn)

# Make executable for job
if process == 0 or process == 3:
  lines = []
  lines.append("#!/bin/sh\n")
  lines.append("export X509_USER_PROXY=$1\n")
  lines.append("voms-proxy-info -all\n")
  lines.append("voms-proxy-info -all -file $1\n")
  lines.append("j=${4:-0}\n")
  lines.append("i=$(($2+5000*$j))\n")
  lines.append("\n")
  lines.append("existname=/eos/user/s/siluo/${5}/out_${i}.root\n")
  lines.append("if [ -e $existname ]\n")
  lines.append("then\n")
  lines.append("  echo 'File already processed'\n")
  lines.append("  exit 0\n")
  lines.append("fi\n")
  lines.append("\n")
  lines.append("# dir=$(pwd)\n")
  lines.append("cd /afs/cern.ch/work/s/siluo/CMSSW/" + cmssw_ver + "/src\n")
  lines.append("export SCRAM_ARCH=slc7_amd64_gcc900\n")
  lines.append("eval `scramv1 runtime -sh`\n")
  lines.append("cd "+curdir+"/\n")
  lines.append("cmsRun NtuplizeSIM.py ifile=$i dataset=$3 outputtag=$5\n")
  lines.append("\n")
  lines.append("filename=/eos/user/s/siluo/${5}/InProgress/out_${i}.root\n")
  lines.append("filesize=$(stat -c %s $filename)\n")
  lines.append("if [ $filesize -gt 1000 ]\n")
  lines.append("then\n")
  lines.append("  mv $filename $existname\n")
  lines.append("  echo 'Job is done...'\n")
  lines.append("else\n")
  lines.append("  rm $filename\n")
  lines.append('  echo "Job failed...Output filesize $filesize. Now retrying......"\n')
  lines.append("  exit 1\n")
  lines.append("fi\n")

  fn = "Submits/MultiSub.sh"
  f = open(fn,"w")
  f.writelines(lines)
  mode = "755"
  os.chmod(fn,int(mode,base=8));
  print('Executable shell script created: '+ fn)

if process == 0 or process == 4:
  lines = []
  lines.append("#!/bin/sh\n")
  lines.append('voms-proxy-init --rfc --voms cms -valid 192:00 -out ${HOME}/x509up\n')
  lines.append('export X509_USER_PROXY=${HOME}/x509up\n')
  lines.append("for i in $(ls *.sub)\n")
  lines.append("do\n")
  lines.append("  echo Submitting $i\n")
  lines.append("  condor_submit $i -batch-name $i\n")
  lines.append("done")

  fn = "Submits/Submission.sh"
  f = open(fn,"w")
  f.writelines(lines)
  mode = "755"
  os.chmod(fn,int(mode,base=8));
  print('Submitting shell script created: '+ fn)

if process == 0 or process == 5:
  lines = []
  lines.append("cmsRun NtuplizeSIM.py ifile=0 dataset={} outputtag={}\n".format(InPath+datasets[0],OutPath+datasets[0]))
  fn = "TestRun.sh"
  f = open(fn,"w")
  f.writelines(lines)
  mode = "755"
  os.chmod(fn,int(mode,base=8));
  print('TestRun shell script created: '+ fn)
