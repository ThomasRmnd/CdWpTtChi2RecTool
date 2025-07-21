import argparse
import sys

parser = argparse.ArgumentParser()

parser.add_argument("--input", type=str, help="Input filepath")
parser.add_argument("--output", type=str, help="Output filepath")
parser.add_argument("--use-simulation", action="store_true", help="Enable SimEvent")

parser.add_argument("--use-auto-factory", action="store_true", help="Use the automatic factory")
parser.add_argument("--water-phase", action="store_true", help="Use water phase algorithms")
parser.add_argument("--manual-reconstruction-mode", type=int, default=1, help="Set the reconstruction mode for the manual factory")
parser.add_argument("--config-file", type=str, default="", help="Path for the config file")

args = parser.parse_args()

ipath = args.input
opath = args.output
use_sim = args.use_simulation

sim_hdr = [
    "/Event/Sim"
]
rec_hdr = [
    "/Event/CdTrackRec", 
    "/Event/WpRec", 
    "/Event/TtRec"
]

# === Sniper ====
import Sniper
Sniper.setLogLevel(1)
task = Sniper.TopTask("task")
task.setLogLevel(1)

# === Profiling ===
import SniperProfiling
prof = task.createSvc("SniperProfiling")
prof.setLogLevel(1)

# === BufferMemMgr ===
import BufferMemMgr
buf_mgr = task.createSvc("BufferMemMgr")
buf_mgr.property("TimeWindow").set([0, 0])

# === Geometry === 
import Geometry
geom = task.createSvc("RecGeomSvc")
geom.property("GeomFile").set("default")
geom.property("GeomPathInRoot").set("JunoGeom")
geom.property("FastInit").set(True)

# === PMTParamSvc ===
pmt_param_svc = task.createSvc("PMTParamSvc")

# === TTGeomSvc ===
tt_geom_svc = task.createSvc("TTGeomSvc")

# === RootIOSvc ===
import RootIOSvc
ifiles = [ipath]
ri_svc = task.createSvc("RootInputSvc/InputSvc")
ri_svc.property("InputFile").set(ifiles)

ofiles = {hdr: opath for hdr in rec_hdr}
if use_sim:
    ofiles.update({hdr: opath for hdr in sim_hdr})

ro_svc = task.createSvc("RootOutputSvc/OutputSvc")
ro_svc.property("OutputStreams").set(ofiles)

# === RecMuonAlg and CdWpTtChi2RecTool ===
import RecMuonAlg
import CdWpTtChi2RecTool
rec_alg = RecMuonAlg.createAlg(task)
rec_alg.setLogLevel(1)
rec_alg.useRecTool("CdWpTtChi2RecTool")

rec_alg.property("Pmt20inchTimeReso").set(8.0)
rec_alg.property("Pmt3inchTimeReso").set(15.0)
rec_alg.property("PmtTTTimeReso").set(2.0)
rec_alg.property("Use3inchPMT").set(True)
rec_alg.property("Use20inchPMT").set(True)
rec_alg.property("ChosenDetectors").set(3) # 1: CD, 2: WP, 4: TT
# rec_alg.property("TimeDiffEventThreshold").set(300.0) # ~200.0 for Joint should be ok
rec_alg.rectool.property("UseAutomaticFactory").set(args.use_auto_factory)
rec_alg.rectool.property("WaterPhase").set(args.water_phase)
rec_alg.rectool.property("ManualReconstructionMode").set(args.manual_reconstruction_mode)
rec_alg.rectool.property("ConfigMap").set(args.config_file)

task.setEvtMax(-1)
# task.show()
if (task.run()):
    print(f"Task finished successfully!")
    sys.exit(0)
else:
    print(f"Task failed!")
    sys.exit(1)