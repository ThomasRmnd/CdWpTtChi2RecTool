import sys

import BufferMemMgr
import Geometry
import RootIOSvc
import Sniper
import SniperProfiling

import RecMuonAlg
import CdWpTtChi2RecTool

ipath = sys.argv[1]
opath = sys.argv[2]

ifile = ipath.split("/")[-1]
ipath = ipath.replace(ifile, "")
ofile = ifile.replace("calib", "rec")

# === Sniper ====
Sniper.setLogLevel(1)
task = Sniper.TopTask("task")
task.setLogLevel(1)

# === Profiling ===
prof = task.createSvc("SniperProfiling")
prof.setLogLevel(1)

# === BufferMemMgr ===
buf_mgr = task.createSvc("BufferMemMgr")
buf_mgr.property("TimeWindow").set([0, 0])

# === Geometry === 
geom = task.createSvc("RecGeomSvc")
geom.property("GeomFile").set("default")
geom.property("GeomPathInRoot").set("JunoGeom")
geom.property("FastInit").set(True)

Sniper.loadDll("libPmtRec.so")
task.property("algs").append("PullSimHeaderAlg")

# === PMTParamSvc ===
pmt_param_svc = task.createSvc("PMTParamSvc")

# === TTGeomSvc ===
tt_geom_svc = task.createSvc("TTGeomSvc")

# === RootIOSvc ===

input_files = [ipath + ifile]

ri_svc = task.createSvc("RootInputSvc/InputSvc")
ri_svc.property("InputFile").set(input_files)

output_files = {
    # === Sim ===
    "/Event/Sim": opath + ofile,
    # === Rec ===
    "/Event/CdTrackRec": opath + ofile,
    "/Event/WpRec": opath + ofile,
    "/Event/TtRec": opath + ofile
}

ro_svc = task.createSvc("RootOutputSvc/OutputSvc")
ro_svc.property("OutputStreams").set(output_files)

# === RecMuonAlg and CdWpTtChi2RecTool ===
import RecMuonAlg
import CdWpTtChi2Rectool
rec_alg = RecMuonAlg.createAlg(task)
rec_alg.setLogLevel(1)
rec_alg.useRecTool("CdWpTtChi2RecTool")

rec_alg.property("ChosenDetectors").set(3) # 1: CD, 2: WP, 4: TT
rec_alg.property("Use20inchPMT").set(True)
rec_alg.property("Use3inchPMT").set(True)
rec_alg.property("Pmt20inchTimeReso").set(8.0)
rec_alg.property("Pmt3inchTimeReso").set(15.0) # 15.0
rec_alg.property("PmtWPTimeReso").set(20.0)
rec_alg.property("PmtTTTimeReso").set(2.0)
rec_alg.property("TimeDiffEventThreshold").set(300.0) # ~200.0 for Joint should be ok

task.setEvtMax(-1)
# task.show()
if (task.run()):
    print("Task finished successfully!")
    sys.exit(0)
else:
    print("Task failed!")
    sys.exit(1)