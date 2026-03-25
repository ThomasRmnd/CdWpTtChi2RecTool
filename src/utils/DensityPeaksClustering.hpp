#ifndef CDWPTTCHI2RECTOOL_UTILS_DENSITYPEAKSCLUSTERING_HPP_
#define CDWPTTCHI2RECTOOL_UTILS_DENSITYPEAKSCLUSTERING_HPP_

#include <cstdlib>
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <vector>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <dirent.h>
#include <sys/stat.h>
#include <numeric>
#include <algorithm>
#include <set>
#include <array>
#include <map>

#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TString.h>
#include <TMath.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TVector3.h>
#include <TStopwatch.h>

#include "RecTools/IRecMuonTool.h"
#include "SniperKernel/AlgFactory.h"
#include "SniperKernel/AlgBase.h"
#include "SniperKernel/ToolBase.h"
#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/SniperDataPtr.h"
#include "RootWriter/RootWriter.h"
#include "Event/WpRecHeader.h"
#include "Geometry/IRecGeomSvc.hh"
#include "Identifier/Identifier.h"
#include "Identifier/JunoDetectorID.h"
#include "Identifier/WpID.h"

#include "utils/RecPmtProp.hpp"

class DensityPeaksClustering {

public:

    void cluster(const RecPmtTable& table);

    const std::vector<double>& getDensities() const;
    const std::vector<double>& getDistances() const;

private:

};

class WpMuonClassifyRecTool :public IRecMuonTool,  public ToolBase {
    
public:
    
    WpMuonClassifyRecTool(const std::string& name);
        
    virtual ~WpMuonClassifyRecTool(); 

    bool reconstruct(RecTrks*);

    void classify(RecPmtTable::const_iterator ftable, RecPmtTable::const_iterator ltable);
    
    bool configure(const Params*, const PmtTable*);

private:
    
    void GroupWPPMT();
    
    void CalcChargePairAve();
    
    void CalcDistance2Cluster();
    
    void InitialEnterPoint();
    
    void InitialExitPoint();
    
    void PrepareCombinationTable();
    
    void generateUniquePairs(const std::vector<int>&, int);
    
    void generateOddUniquePairs(const std::vector<int>&, int);
    
    bool isValidSet(const std::vector<std::pair<int, int>>&, const std::vector<int>&);
    
    double CalcAngle(const std::vector<TVector3>&, int);
    
    bool LoadAdditionalGainCorrection();

    const PmtTable* m_ptable;
    int TotalWPPMT;
    std::vector<TVector3> ALL_WPPMT_pos;
    std::vector<std::vector<int> > ALL_WPPMT_pair;
    std::vector<int> ALL_WPPMT_npair;
    std::vector<int> ALL_WPPMT_highchargepmt_in_group;
    std::vector<double> ALL_WPPMT_Charge;
    std::vector<double> ALL_WPPMT_HitTime;
    std::vector<double> ALL_WPPMT_AveCharge;
    std::vector<double> ALL_WPPMT_Distance;
    bool bUseAdditionalGainCorrection = false;
    std::string AdditionalGainCorrectionPath;
    std::vector<double> ALL_WPPMT_gaincorrection;
    double WPmaxpulsecharge = 0.;
    double totalWPcharge = 0.;
    std::vector<int> ClusterCenterPMTID;
    std::vector<int> ClusterHighestChargePMTID;
    std::vector<int> EnterPMTID;
    std::vector<int> ExitPMTID;
    static const int nlongdis = 5;
    std::vector<double> longdisvec = std::vector<double>(nlongdis, 0.);
    double longdisarray[nlongdis];
    double initialchargecut = 28.;
    double WPmaxpulsechargeth = 200;
    double disth = 6500.;
    const double WPPMT_R = 20500.;
    int recoclassid;
    int nrecocluster;
    int initialenterpmtid;
    double initialenterx;
    double initialentery;
    double initialenterz;
    int initialexitpmtid;
    double initialexitx;
    double initialexity;
    double initialexitz;
    static const int nmaxenter = 20;
    static const int nmaxexit = 20;
    int nenter;
    int nexit;
    double recoenterx[nmaxenter];
    double recoentery[nmaxenter];
    double recoenterz[nmaxenter];
    double recoexitx[nmaxexit];
    double recoexity[nmaxexit];
    double recoexitz[nmaxexit];
    static const int maxuniquesets = 10;
    std::vector<std::vector<std::vector<std::pair<int, int> > > > uniqueSetsvec;
        
};

#endif // CDWPTTCHI2RECTOOL_UTILS_DENSITYPEAKSCLUSTERING_HPP_