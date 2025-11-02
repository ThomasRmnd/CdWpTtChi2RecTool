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

parser.add_argument("--time-window", nargs=2, type=float, metavar=("START", "END"), default=(-1e-6, 1e-6), help="Buffer time window")
parser.add_argument("--log-level", type=int, default=3, help="Log level (default: 1)")

args = parser.parse_args()

ipath = args.input
opath = args.output

lower_tw, upper_tw = args.time_window
loglevel = args.log_level

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
Sniper.setLogLevel(loglevel)
task = Sniper.TopTask("task")
task.setLogLevel(loglevel)

# === Profiling ===
import SniperProfiling
prof = task.createSvc("SniperProfiling")
prof.setLogLevel(loglevel)

# === BufferMemMgr ===
import BufferMemMgr
buf_mgr = task.createSvc("BufferMemMgr")
buf_mgr.property("TimeWindow").set([lower_tw, upper_tw])

# === Geometry === 
import Geometry
geom = task.createSvc("RecGeomSvc")
geom.property("GeomFile").set("default")
geom.property("GeomPathInRoot").set("JunoGeom")
geom.property("FastInit").set(True)
pmt_svc = task.createSvc("PMTParamSvc")
ttg_svc = task.createSvc("TTGeomSvc")

# === RootIOSvc ===
import RootIOSvc
ri_svc = task.createSvc("RootInputSvc/InputSvc")
ri_svc.property("InputFile").set(ipath)

ofiles = {hdr: opath for hdr in rec_hdr}
if args.use_simulation:
    ofiles.update({hdr: opath for hdr in sim_hdr})

ro_svc = task.createSvc("RootOutputSvc/OutputSvc")
ro_svc.property("OutputStreams").set(ofiles)

# === RecMuonAlg and CdWpTtChi2RecTool ===
import RecMuonAlg
import CdWpTtChi2RecTool
alg = RecMuonAlg.createAlg(task)
alg.setLogLevel(loglevel)

alg.useLoader("JointLoader")
alg.loader.property("TimeWindow").set([-500.0, 500.0]) # ns

alg.useCdFiller("CdRangeFiller")
alg.cdfiller.property("Pmt3inchTimeReso").set(15.0)
alg.cdfiller.property("Pmt20inchTimeReso").set(8.0)

alg.useWpFiller("WpRangeFiller")
alg.wpfiller.property("PmtTimeReso").set(8.0)

alg.useTtFiller("TtRangeFiller")
alg.ttfiller.property("PmtTimeReso").set(2.0)

alg.useRecTool("CdWpTtChi2RecTool")
alg.rectool.property("UseAutomaticFactory").set(args.use_auto_factory)
alg.rectool.property("WaterPhase").set(args.water_phase)
alg.rectool.property("ManualReconstructionMode").set(args.manual_reconstruction_mode)
alg.rectool.property("ConfigFile").set(args.config_file)

task.setEvtMax(-1)
if (task.run()):
    print(f"Task finished successfully!")
    sys.exit(0)
else:
    print(f"Task failed!")
    sys.exit(1)