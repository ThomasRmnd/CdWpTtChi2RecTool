#include "utils/DensityPeaksClustering.hpp"

WpMuonClassifyRecTool::WpMuonClassifyRecTool(const std::string& name) : 
    ToolBase(name), 
    m_ptable(NULL)
{
    declProp("WpMuonClassifyRecToolInitialChargeCut", initialchargecut = 28.);
    declProp("WpMuonClassifyRecToolMaxChargeThreshold", WPmaxpulsechargeth = 200.);
    declProp("WpMuonClassifyRecToolDistanceThreshold", disth = 6500.);
    declProp("UseAdditionalGainCorrection", bUseAdditionalGainCorrection = false);
    declProp("AdditionalGainCorrectionPath", AdditionalGainCorrectionPath);
}

bool WpMuonClassifyRecTool::configure(const Params* pars, const PmtTable* ptab) {
    m_ptable = ptab;

    LogInfo << "TotalWPPMT = " << TotalWPPMT << ", initialchargecut = " << initialchargecut << ", WPmaxpulsechargeth = " << WPmaxpulsechargeth << ", disth = " << disth << std::endl;
    GroupWPPMT();
    PrepareCombinationTable();
    ALL_WPPMT_gaincorrection.resize(TotalWPPMT, 1.);
    if (bUseAdditionalGainCorrection) {
        bool bload = LoadAdditionalGainCorrection();
        if (!bload) return false;
    }
    return true; 
}
bool WpMuonClassifyRecTool::LoadAdditionalGainCorrection()
{
    FILE *finputtxt = NULL;
    int copyno;
    double gainratio;
    finputtxt = fopen(Form("%s", AdditionalGainCorrectionPath.c_str()), "r");
    if (finputtxt==NULL) {
        LogError << AdditionalGainCorrectionPath << " is not found." << std::endl;
        return false;
    } else {
        LogInfo << "Reading: " << AdditionalGainCorrectionPath << std::endl;
    }
    while(fscanf(finputtxt, "%d\t%lf", &copyno, &gainratio)!=EOF) {
        if (copyno>=kOFFSET_WP_PMT) copyno = copyno - kOFFSET_WP_PMT;
        ALL_WPPMT_gaincorrection.at(copyno) = gainratio;
        if (copyno<10) LogInfo << "copyno = " << copyno << ", gainratio = " << gainratio<< std::endl;
    }
    fclose(finputtxt);
    return true;
}

void WpMuonClassifyRecTool::classify(RecPmtTable::const_iterator ftable, RecPmtTable::const_iterator ltable) {
    
    // GroupWPPMT ==> create vector of neighboring PMTs with radius 2500 mm
    double distcut = 2500.0;
    std::vector<RecPmtTable> table_pairs(std::distance(ftable, ltable));
    for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
        for (RecPmtTable::const_iterator jt = ftable; jt != ltable; ++jt) {
            if (it == jt) continue;
            double dist = mag(jt->pos - it->pos);
            if (distcut < dist) continue;
            table_pairs[std::distance(ftable, it)].push_back(*jt);
        }
    }

    // Count the number of time a PMT is neighboring a high charge PMT
    std::vector<int> table_nb_neighboring_highq(std::distance(ftable, ltable), 0);
    for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
        if (it->totq < initialchargecut) continue;
        RecPmtTable::const_iterator jftable = table_pairs[std::distance(ftable, it)].begin();
        RecPmtTable::const_iterator jltable = table_pairs[std::distance(ftable, it)].end();
        for (RecPmtTable::const_iterator jt = jftable; jt != jltable; ++jt) {
            table_nb_neighboring_highq[std::distance(jftable, jt)]++;
        }
    }

    // Get the highest charged PMT
    RecPmtTable::const_iterator it_highq = ltable;
    double maxq = 0.0;
    for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
        if (it->totq <= maxq) continue;
        maxq = it->totq;
        it_highq = it;
    }

    std::vector<double> table_averageq(std::distance(ftable, ltable), 0.0);
    double reductionq = 10.0;
    for (RecPmtTable::const_iterator it = ftable; it != ltable; ++it) {
        table_averageq[std::distance(ftable, it)] = it->totq / table_nb_neighboring_highq[std::distance(ftable, it)];
        RecPmtTable::const_iterator jftable = table_pairs[std::distance(ftable, it)].begin();
        RecPmtTable::const_iterator jltable = table_pairs[std::distance(ftable, it)].end();
        for (RecPmtTable::const_iterator jt = jftable; jt != jltable; ++jt) {
            table_averageq[std::distance(jftable, jt)] += jt->totq / table_nb_neighboring_highq[std::distance(ftable, it)];
        }
        if (table_averageq[std::distance(ftable, it)] >= it->totq / reductionq) continue;
        table_averageq[std::distance(ftable, it)] = it->totq / reductionq;
    }
}

bool WpMuonClassifyRecTool::reconstruct(RecTrks* trks) {
    const PmtTable& ptab = *m_ptable;
    WPmaxpulsecharge = 0.;
    totalWPcharge = 0.;
    int WPmaxpulsechargepmtid = -1;
    int track_flag = -1;
    ALL_WPPMT_Charge = std::vector<double>(TotalWPPMT);
    ALL_WPPMT_HitTime = std::vector<double>(TotalWPPMT);
    ALL_WPPMT_AveCharge = std::vector<double>(TotalWPPMT);
    ALL_WPPMT_Distance = std::vector<double>(TotalWPPMT);
    ALL_WPPMT_highchargepmt_in_group = std::vector<int>(TotalWPPMT, 0);

    for (int i = 0; i < int(ptab.size()); i++) {
        ALL_WPPMT_Charge.at(ptab.at(i).pmtid) = ptab.at(i).q*ALL_WPPMT_gaincorrection.at(ptab.at(i).pmtid);
        ALL_WPPMT_HitTime.at(ptab.at(i).pmtid) = ptab.at(i).fht;
        if ((ptab.at(i).q*ALL_WPPMT_gaincorrection.at(ptab.at(i).pmtid))>initialchargecut) {
            for (int j = 0; j < ALL_WPPMT_pair.at(ptab.at(i).pmtid).size(); j++) {
                ALL_WPPMT_highchargepmt_in_group.at(ALL_WPPMT_pair.at(ptab.at(i).pmtid).at(j))++;
            }
        }
        if (WPmaxpulsecharge<ptab.at(i).q*ALL_WPPMT_gaincorrection.at(ptab.at(i).pmtid)) {
            WPmaxpulsecharge = ptab.at(i).q*ALL_WPPMT_gaincorrection.at(ptab.at(i).pmtid);
            WPmaxpulsechargepmtid = ptab.at(i).pmtid;
        }
    }

    CalcChargePairAve();
    CalcDistance2Cluster();
    if (WPmaxpulsecharge<WPmaxpulsechargeth) {
        recoclassid = 0; 
        track_flag = -3;
    } else if (nrecocluster==0) {
        recoclassid = 0; 
        track_flag = -2;
    } else if (nrecocluster==1) {
        recoclassid = 0; 
        track_flag = -1;
    } else if (nrecocluster==2) {
        recoclassid = 1; 
    } else {
        recoclassid = 2; 
    }
    TVector3 start, end;
    std::vector<TVector3> multistart;
    std::vector<TVector3> multiend;
    std::vector<TVector3> multitrack;
    int ntracks;
    double minangle = 999999.;
    double tmpangle;
    int minangleii = 0;
    if (recoclassid==1) {
        InitialEnterPoint();
        InitialExitPoint();
        start.SetXYZ(initialenterx, initialentery, initialenterz);
        end.SetXYZ(initialexitx, initialexity, initialexitz);
        if ((ALL_WPPMT_pos.at(ClusterCenterPMTID.at(0)).Z())>(ALL_WPPMT_pos.at(ClusterCenterPMTID.at(1)).Z())) {
            trks->addTrk(start, end, ALL_WPPMT_HitTime.at(ClusterCenterPMTID.at(0)), ALL_WPPMT_HitTime.at(ClusterCenterPMTID.at(1)), totalWPcharge, double(0), 0);
        } else {
            trks->addTrk(start, end, ALL_WPPMT_HitTime.at(ClusterCenterPMTID.at(1)), ALL_WPPMT_HitTime.at(ClusterCenterPMTID.at(0)), totalWPcharge, double(0), 0);
        }
        nenter = 1;
        recoenterx[0] = initialenterx; 
        recoentery[0] = initialentery; 
        recoenterz[0] = initialenterz; 
        nexit = 1;
        recoexitx[0] = initialexitx; 
        recoexity[0] = initialexity; 
        recoexitz[0] = initialexitz; 
    } else if (recoclassid==0) {
        initialenterx = ALL_WPPMT_pos.at(WPmaxpulsechargepmtid).X();
        initialentery = ALL_WPPMT_pos.at(WPmaxpulsechargepmtid).Y();
        initialenterz = ALL_WPPMT_pos.at(WPmaxpulsechargepmtid).Z();
        start.SetXYZ(initialenterx, initialentery, initialenterz);
        end.SetXYZ(-999999., -999999., -999999.);
        trks->addTrk(start, end, ALL_WPPMT_HitTime.at(WPmaxpulsechargepmtid), -999999., totalWPcharge, double(track_flag), track_flag);
        nenter = 1;
        recoenterx[0] = initialenterx; 
        recoentery[0] = initialentery; 
        recoenterz[0] = initialenterz; 
        nexit = 0;

    } else if (recoclassid==2) {
        if ((nrecocluster)%2==1) {
            multistart.resize((nrecocluster+1)/2);
            multiend.resize((nrecocluster+1)/2);
            ntracks = ((nrecocluster+1)/2);
            multitrack.resize((nrecocluster+1)/2);
        } else {
            multistart.resize((nrecocluster)/2);
            multiend.resize((nrecocluster)/2);
            ntracks = ((nrecocluster)/2);
            multitrack.resize((nrecocluster)/2);
        }
        int tmppmtid;
        for (int i = 0; i < int(ClusterCenterPMTID.size()); i++) {
            for (int j = i+1; j < int(ClusterCenterPMTID.size()); j++) {
                if (ALL_WPPMT_pos.at(ClusterCenterPMTID.at(i)).Z()<ALL_WPPMT_pos.at(ClusterCenterPMTID.at(j)).Z()) {
                    tmppmtid = ClusterCenterPMTID.at(i);
                    ClusterCenterPMTID.at(i) = ClusterCenterPMTID.at(j);
                    ClusterCenterPMTID.at(j) = tmppmtid;
                }
            }
        }
        if (nrecocluster<=8) {
            nenter = ntracks;
            nexit = ntracks;
            for (int ii = 0; ii < uniqueSetsvec.at(nrecocluster).size(); ii++) {
                for (int jj = 0; jj < ntracks; jj++) {
                    multistart.at(jj) = ALL_WPPMT_pos.at(ClusterCenterPMTID.at(uniqueSetsvec.at(nrecocluster).at(ii).at(jj).first));
                    multiend.at(jj) = ALL_WPPMT_pos.at(ClusterCenterPMTID.at(uniqueSetsvec.at(nrecocluster).at(ii).at(jj).second));
                    multitrack.at(jj) = multiend.at(jj)-multistart.at(jj);
                }
                tmpangle = CalcAngle(multitrack, ntracks);
                if (minangle>tmpangle) {
                    minangle=tmpangle;
                    minangleii = ii;
                }
            }
            for (int jj = 0; jj < ntracks; jj++) {
                trks->addTrk(ALL_WPPMT_pos.at(ClusterCenterPMTID.at(uniqueSetsvec.at(nrecocluster).at(minangleii).at(jj).first)), ALL_WPPMT_pos.at(ClusterCenterPMTID.at(uniqueSetsvec.at(nrecocluster).at(minangleii).at(jj).second)), ALL_WPPMT_HitTime.at(ClusterCenterPMTID.at(uniqueSetsvec.at(nrecocluster).at(minangleii).at(jj).first)), ALL_WPPMT_HitTime.at(ClusterCenterPMTID.at(uniqueSetsvec.at(nrecocluster).at(minangleii).at(jj).second)), totalWPcharge, double(jj), jj);
                recoenterx[jj] = ALL_WPPMT_pos.at(ClusterCenterPMTID.at(uniqueSetsvec.at(nrecocluster).at(minangleii).at(jj).first)).X();
                recoentery[jj] = ALL_WPPMT_pos.at(ClusterCenterPMTID.at(uniqueSetsvec.at(nrecocluster).at(minangleii).at(jj).first)).Y();
                recoenterz[jj] = ALL_WPPMT_pos.at(ClusterCenterPMTID.at(uniqueSetsvec.at(nrecocluster).at(minangleii).at(jj).first)).Z();
                recoexitx[jj] = ALL_WPPMT_pos.at(ClusterCenterPMTID.at(uniqueSetsvec.at(nrecocluster).at(minangleii).at(jj).second)).X();
                recoexity[jj] = ALL_WPPMT_pos.at(ClusterCenterPMTID.at(uniqueSetsvec.at(nrecocluster).at(minangleii).at(jj).second)).Y();
                recoexitz[jj] = ALL_WPPMT_pos.at(ClusterCenterPMTID.at(uniqueSetsvec.at(nrecocluster).at(minangleii).at(jj).second)).Z();
            }
        } else {
            if (ClusterCenterPMTID.size()%2==1) {
                ClusterCenterPMTID.push_back(ClusterCenterPMTID.at(ClusterCenterPMTID.size()-1));
            }
            nenter = ClusterCenterPMTID.size()/2;
            nexit = ClusterCenterPMTID.size()/2;
            for (int jj = 0; jj < ClusterCenterPMTID.size()/2; jj++) {
                trks->addTrk(ALL_WPPMT_pos.at(ClusterCenterPMTID.at(jj)), ALL_WPPMT_pos.at(ClusterCenterPMTID.at(jj+nenter)), ALL_WPPMT_HitTime.at(ClusterCenterPMTID.at(jj)), ALL_WPPMT_HitTime.at(ClusterCenterPMTID.at(jj+nenter)), totalWPcharge, double(jj), jj);
                recoenterx[jj] = ALL_WPPMT_pos.at(ClusterCenterPMTID.at(jj)).X();
                recoentery[jj] = ALL_WPPMT_pos.at(ClusterCenterPMTID.at(jj)).Y();
                recoenterz[jj] = ALL_WPPMT_pos.at(ClusterCenterPMTID.at(jj)).Z();
                recoexitx[jj] = ALL_WPPMT_pos.at(ClusterCenterPMTID.at(jj+nenter)).X();
                recoexity[jj] = ALL_WPPMT_pos.at(ClusterCenterPMTID.at(jj+nenter)).Y();
                recoexitz[jj] = ALL_WPPMT_pos.at(ClusterCenterPMTID.at(jj+nenter)).Z();
            }
        }
    } else {
        initialenterx = -999999;
        initialentery = -999999;
        initialenterz = -999999;
        initialexitx = -999999;
        initialexity = -999999;
        initialexitz = -999999;
    }

    return true;
}
void WpMuonClassifyRecTool::GroupWPPMT()
{
    double dis; 
    const double discut = 2500.;
    ALL_WPPMT_pair.resize(TotalWPPMT);
    ALL_WPPMT_npair.resize(TotalWPPMT, 0);
    double mindisperpmt = 999999.;
    for (int i = 0; i < TotalWPPMT; i++) {
        mindisperpmt = 999999.;
        for (int j = 0; j < TotalWPPMT; j++) {
            if (i==j) continue;
            dis = sqrt((ALL_WPPMT_pos.at(i).X()-ALL_WPPMT_pos.at(j).X())*(ALL_WPPMT_pos.at(i).X()-ALL_WPPMT_pos.at(j).X())+(ALL_WPPMT_pos.at(i).Y()-ALL_WPPMT_pos.at(j).Y())*(ALL_WPPMT_pos.at(i).Y()-ALL_WPPMT_pos.at(j).Y())+(ALL_WPPMT_pos.at(i).Z()-ALL_WPPMT_pos.at(j).Z())*(ALL_WPPMT_pos.at(i).Z()-ALL_WPPMT_pos.at(j).Z()));
            if (dis<discut) {
                ALL_WPPMT_pair.at(i).push_back(j);
                ALL_WPPMT_npair.at(i)++;
            }    
            if ((dis<mindisperpmt)&&(i!=j)) {
                mindisperpmt = dis; 
            }    
        }   
    }
}
void WpMuonClassifyRecTool::CalcChargePairAve()
{
    const double reductionth = 10.;
    for (int i = 0; i < TotalWPPMT; i++) {
        if (ALL_WPPMT_Charge.at(i)==0) continue;
        ALL_WPPMT_AveCharge.at(i) = ALL_WPPMT_Charge.at(i)/ALL_WPPMT_npair.at(i);
        for (int j = 0; j < ALL_WPPMT_npair.at(i); j++) {
            ALL_WPPMT_AveCharge.at(i) += ALL_WPPMT_Charge.at(ALL_WPPMT_pair.at(i).at(j))/ALL_WPPMT_npair.at(i);
        }
        if (ALL_WPPMT_AveCharge.at(i)<(ALL_WPPMT_Charge.at(i)/reductionth)) {
            ALL_WPPMT_AveCharge.at(i) = ALL_WPPMT_Charge.at(i)/reductionth;
        }
    }
    return;
}
void WpMuonClassifyRecTool::CalcDistance2Cluster() 
{
    ALL_WPPMT_Distance = std::vector<double>(TotalWPPMT, -1.);
    double tmpmindis;
    double tmpmaxdis = 0; 
    double dis; 
    double maxchargeave = -1.; 
    int maxchargepmtid = -1;
    TVector3 tmpvec;
    int disrank;
    int pairedpmtid;
    nrecocluster = 0; 
    longdisvec = std::vector<double>(nlongdis, 0.); 
    ClusterCenterPMTID = std::vector<int>();
    ClusterHighestChargePMTID = std::vector<int>();
    double highestchargecluster{0};
    double highestchargeclusterpmtid{0};
    std::vector<int> misspmtid = std::vector<int>();

    for (int i = 0; i < TotalWPPMT; i++) {
        tmpmindis = 999999.;
        if (maxchargeave<ALL_WPPMT_AveCharge.at(i)) {
            maxchargeave = ALL_WPPMT_AveCharge.at(i);
            maxchargepmtid = i; 
        } else if (maxchargeave==ALL_WPPMT_AveCharge.at(i)) {
            if (ALL_WPPMT_Charge.at(maxchargepmtid)<ALL_WPPMT_Charge.at(i)) {
                maxchargeave = ALL_WPPMT_AveCharge.at(i);
                maxchargepmtid = i; 
            }    
        }    
        if (ALL_WPPMT_Charge.at(i)<initialchargecut) continue;
        if ((ALL_WPPMT_highchargepmt_in_group.at(i)<2)||(ALL_WPPMT_highchargepmt_in_group.at(i)<int(0.34*ALL_WPPMT_npair.at(i)))) continue;
        for (int j = 0; j < TotalWPPMT; j++) {
            if (i==j) continue;
            if (ALL_WPPMT_Charge.at(j)<initialchargecut) continue;
            if ((ALL_WPPMT_highchargepmt_in_group.at(j)<2)||(ALL_WPPMT_highchargepmt_in_group.at(j)<int(0.34*ALL_WPPMT_npair.at(j)))) continue;
            if (ALL_WPPMT_AveCharge.at(i)<ALL_WPPMT_AveCharge.at(j)) {
                tmpvec = ALL_WPPMT_pos.at(j);
                tmpvec.RotateZ(-ALL_WPPMT_pos.at(i).Phi());
                tmpvec.RotateY(-ALL_WPPMT_pos.at(i).Theta());
                dis = abs(WPPMT_R*tmpvec.Theta());
                if (dis<tmpmindis) {
                    pairedpmtid = j; 
                    tmpmindis = dis; 
                }    
            } else if (ALL_WPPMT_AveCharge.at(i)==ALL_WPPMT_AveCharge.at(j)) {
                if (ALL_WPPMT_Charge.at(i)<ALL_WPPMT_Charge.at(j)) {
                    tmpvec = ALL_WPPMT_pos.at(j);
                    tmpvec.RotateZ(-ALL_WPPMT_pos.at(i).Phi());
                    tmpvec.RotateY(-ALL_WPPMT_pos.at(i).Theta());
                    dis = abs(WPPMT_R*tmpvec.Theta());
                    if (dis<tmpmindis) {
                        pairedpmtid = j;
                        tmpmindis = dis;
                    }
                }
            }
        }
        ALL_WPPMT_Distance.at(i) = tmpmindis;
        if (ALL_WPPMT_Distance.at(i)==999999.) misspmtid.push_back(i);
    }
    for (int i = 0; i < TotalWPPMT; i++) {
        if (ALL_WPPMT_Charge.at(i)<initialchargecut) continue;
        if ((ALL_WPPMT_highchargepmt_in_group.at(i)<2)||(ALL_WPPMT_highchargepmt_in_group.at(i)<int(0.34*ALL_WPPMT_npair.at(i)))) continue;
        tmpvec = ALL_WPPMT_pos.at(i);
        tmpvec.RotateZ(-ALL_WPPMT_pos.at(maxchargepmtid).Phi());
        tmpvec.RotateY(-ALL_WPPMT_pos.at(maxchargepmtid).Theta());
        dis = abs(WPPMT_R*tmpvec.Theta());
        if (dis>tmpmaxdis) {
            tmpmaxdis = dis;
        }
    }
    ALL_WPPMT_Distance.at(maxchargepmtid) = tmpmaxdis;
    for (int i = 0; i < int(misspmtid.size()); i++) {
        if (misspmtid.at(i)==maxchargepmtid) continue;
        tmpmindis = 999999.;
        for (int j = 0; j < TotalWPPMT; j++) {
            if (misspmtid.at(i)==j) continue;
            if (ALL_WPPMT_AveCharge.at(misspmtid.at(i))<ALL_WPPMT_AveCharge.at(j)) {
                tmpvec = ALL_WPPMT_pos.at(j);
                tmpvec.RotateZ(-ALL_WPPMT_pos.at(misspmtid.at(i)).Phi());
                tmpvec.RotateY(-ALL_WPPMT_pos.at(misspmtid.at(i)).Theta());
                dis = abs(WPPMT_R*tmpvec.Theta());
                if (dis<tmpmindis) {
                    pairedpmtid = j;
                    tmpmindis = dis;
                }
            } else if (ALL_WPPMT_AveCharge.at(misspmtid.at(i))==ALL_WPPMT_AveCharge.at(j)) {
                if (ALL_WPPMT_Charge.at(misspmtid.at(i))<ALL_WPPMT_Charge.at(j)) {

                    tmpvec = ALL_WPPMT_pos.at(j);
                    tmpvec.RotateZ(-ALL_WPPMT_pos.at(misspmtid.at(i)).Phi());
                    tmpvec.RotateY(-ALL_WPPMT_pos.at(misspmtid.at(i)).Theta());
                    dis = abs(WPPMT_R*tmpvec.Theta());
                    if (dis<tmpmindis) {
                        pairedpmtid = j;
                        tmpmindis = dis;
                    }
                }
            }
        }
        ALL_WPPMT_Distance.at(misspmtid.at(i)) = tmpmindis;
    }
    for (int i = 0; i < TotalWPPMT; i++) {
        disrank = -1;
        for (int j = 0; j < nlongdis; j++) {
            if (longdisvec.at(j)<ALL_WPPMT_Distance.at(i)) {
                disrank = j;
                break;
            }
        }
        if (disrank>=0) {
            for (int j = nlongdis-1; j > 0; j--) {
                if (j==disrank) break;
                longdisvec.at(j) = longdisvec.at(j-1);
            }
            longdisvec.at(disrank) = ALL_WPPMT_Distance.at(i);
        }
        if (ALL_WPPMT_Distance.at(i)>disth) {
            ClusterCenterPMTID.push_back(i);
            nrecocluster++;
            highestchargecluster = 0.;
            for (int j = 0; j < ALL_WPPMT_npair.at(i); j++) {
                if (highestchargecluster<ALL_WPPMT_Charge.at(ALL_WPPMT_pair.at(i).at(j))) {
                    highestchargecluster = ALL_WPPMT_Charge.at(ALL_WPPMT_pair.at(i).at(j));
                    highestchargeclusterpmtid = ALL_WPPMT_pair.at(i).at(j);
                }
            }
            ClusterHighestChargePMTID.push_back(highestchargeclusterpmtid);
        }
    }
    for (int j = 0; j < nlongdis; j++) {
        longdisarray[j] = longdisvec.at(j);
    }
    return;
}
void WpMuonClassifyRecTool::InitialEnterPoint()
{
    int highestzpmtid{0};
    double highestz = -999999.;
    TVector3 tmpvec;
    TVector3 tmpinitialenter;
    tmpinitialenter.SetXYZ(0, 0, 0);
    double dis{0}; 
    int initialenterpmtid{0};
    for (int i = 0; i < int(ClusterCenterPMTID.size()); i++) {
        if (ALL_WPPMT_pos.at(ClusterCenterPMTID.at(i)).Z()>highestz) {
            highestz = ALL_WPPMT_pos.at(ClusterCenterPMTID.at(i)).Z();
            highestzpmtid = ClusterCenterPMTID.at(i);
        }    
    }    
    double highestcharge = 0; 
    for (int i = 0; i < TotalWPPMT; i++) {
        if (ALL_WPPMT_Charge.at(i)==0) continue;
        tmpvec = ALL_WPPMT_pos.at(i);
        tmpvec.RotateZ(-ALL_WPPMT_pos.at(highestzpmtid).Phi());
        tmpvec.RotateY(-ALL_WPPMT_pos.at(highestzpmtid).Theta());

        dis = abs(WPPMT_R*tmpvec.Theta());
        if ((dis<disth)&&(highestcharge<ALL_WPPMT_Charge.at(i))) {
            initialenterpmtid = i; 
            highestcharge = ALL_WPPMT_Charge.at(i);
        }    
        if (dis<disth) {
            tmpinitialenter += (ALL_WPPMT_Charge.at(i)*ALL_WPPMT_Charge.at(i))*ALL_WPPMT_pos.at(i);
        }    
    } 
    tmpinitialenter = WPPMT_R/tmpinitialenter.Mag()*tmpinitialenter;
    initialenterx = tmpinitialenter.X();
    initialentery = tmpinitialenter.Y();
    initialenterz = tmpinitialenter.Z();
}
void WpMuonClassifyRecTool::InitialExitPoint()
{   
    TVector3 tmpvec;
    double dis{0};
    int lowestzpmtid{0};
    double lowestz = 999999.;
    TVector3 tmpinitialexit; 
    tmpinitialexit.SetXYZ(0, 0, 0);
    int initialexitpmtid2{0};
    for (int i = 0; i < int(ClusterCenterPMTID.size()); i++) {
        if (ALL_WPPMT_pos.at(ClusterCenterPMTID.at(i)).Z()<lowestz) {
            lowestz = ALL_WPPMT_pos.at(ClusterCenterPMTID.at(i)).Z();
            lowestzpmtid = ClusterCenterPMTID.at(i);
        }
    }
    double highestcharge = 0;
    for (int i = 0; i < TotalWPPMT; i++) {
        if (ALL_WPPMT_Charge.at(i)==0) continue;
        tmpvec = ALL_WPPMT_pos.at(i);
        tmpvec.RotateZ(-ALL_WPPMT_pos.at(lowestzpmtid).Phi());
        tmpvec.RotateY(-ALL_WPPMT_pos.at(lowestzpmtid).Theta());
        dis = abs(WPPMT_R*tmpvec.Theta());
        if ((dis<disth)&&(highestcharge<ALL_WPPMT_Charge.at(i))&&(ALL_WPPMT_pos.at(i).Z()<initialenterz)) {
            initialexitpmtid2 = i;
            highestcharge = ALL_WPPMT_Charge.at(i);
        }
        if (dis<disth) {
            tmpinitialexit += (ALL_WPPMT_Charge.at(i))*ALL_WPPMT_pos.at(i);
        }
    }
    initialexitx = ALL_WPPMT_pos.at(initialexitpmtid2).X();
    initialexity = ALL_WPPMT_pos.at(initialexitpmtid2).Y();
    initialexitz = ALL_WPPMT_pos.at(initialexitpmtid2).Z();
}
double WpMuonClassifyRecTool::CalcAngle(const std::vector<TVector3>& multitrack, int ntracks)
{
    double angle = 0.;
    if (ntracks==2) {
        angle = abs(multitrack.at(0).Angle(multitrack.at(1)));
    } else if (ntracks==3) {
        angle += abs(multitrack.at(0).Angle(multitrack.at(1)));
        angle += abs(multitrack.at(1).Angle(multitrack.at(2)));
        angle += abs(multitrack.at(0).Angle(multitrack.at(2)));
    } else if (ntracks==4) {
        angle += abs(multitrack.at(0).Angle(multitrack.at(1)));
        angle += abs(multitrack.at(0).Angle(multitrack.at(2)));
        angle += abs(multitrack.at(0).Angle(multitrack.at(3)));
        angle += abs(multitrack.at(1).Angle(multitrack.at(2)));
        angle += abs(multitrack.at(1).Angle(multitrack.at(3)));

        angle += abs(multitrack.at(2).Angle(multitrack.at(3)));
    }
    return angle;
}
void WpMuonClassifyRecTool::PrepareCombinationTable()
{
    uniqueSetsvec.resize(maxuniquesets);
    std::vector<int> balls;
    for (int i = 3; i < 9; i++) {
        if (i%2==0) {
            balls = std::vector<int>();
            for (int j = 0; j < i; j++) {
                balls.push_back(j);
            }
            generateUniquePairs(balls, i);
        } else {
            for (int j = 0; j < i-1; j++) {
                balls = std::vector<int>();
                for (int k = 0; k < i+1; k++) {
                    balls.push_back(k);
                }
                balls.at(balls.size()-1) = j;
                std::sort(balls.begin(), balls.end());
                generateOddUniquePairs(balls, i);
            }
        }
    }
}
// Based on the code from ChatGPT
void WpMuonClassifyRecTool::generateUniquePairs(const std::vector<int>& balls, int num) {
    std::vector<std::pair<int, int>> pairs;
    // Generate all possible pairs from the balls
    for (size_t i = 0; i < balls.size(); ++i) {
        for (size_t j = i + 1; j < balls.size(); ++j) {
            pairs.emplace_back(balls[i], balls[j]);
        }   
    }   
    // To store unique sets of pairs
    std::set<std::set<std::pair<int, int> > > uniqueSets;
    // Generate combinations of 3 pairs
    std::vector<int> comb(pairs.size(), 0); 
    std::fill(comb.end() - balls.size()/2, comb.end(), 1); 
    do {
        std::vector<std::pair<int, int> > selectedPairs;
        std::set<int> usedBalls;
        for (size_t i = 0; i < pairs.size(); ++i) {
            if (comb[i]) {
                selectedPairs.push_back(pairs[i]);
                usedBalls.insert(pairs[i].first);
                usedBalls.insert(pairs[i].second);
            }   
        }   
        // Ensure no ball is used more than once
        if (usedBalls.size() == balls.size()) {
            // Insert the set of pairs into the set of unique combinations
            std::set<std::pair<int, int>> combination(selectedPairs.begin(), selectedPairs.end());
            uniqueSets.insert(combination);
        }   
    } while (std::next_permutation(comb.begin(), comb.end()));
    int ii = 0;
    int jj = 0;
    int npatterns = 0;

    uniqueSetsvec.at(num).resize(uniqueSets.size());
    for (const auto& uset : uniqueSets) {
        uniqueSetsvec.at(num).at(ii).resize(uset.size());
        jj = 0;
        for (const auto& p : uset) {
            uniqueSetsvec.at(num).at(ii).at(jj).first = p.first;
            uniqueSetsvec.at(num).at(ii).at(jj).second = p.second;
            jj++;
        }
        npatterns++;
        ii++;
    }
}
void WpMuonClassifyRecTool::generateOddUniquePairs(const std::vector<int>& balls, int num)
{
    // Step 1: Define the balls
    //    vector<int> balls = {0, 0, 1, 2, 3, 4};
    // Step 2: Generate all pairs, excluding pairs with the same value
    std::vector<std::pair<int, int>> all_pairs;
    for (size_t i = 0; i < balls.size(); ++i) {
        for (size_t j = i + 1; j < balls.size(); ++j) {
            if (balls[i] != balls[j]) { // Exclude pairs like (1, 1)
                all_pairs.emplace_back(balls[i], balls[j]);
            }   
        }   
    }   
    // Step 3: Generate all unique sets of 3 pairs
    std::set<std::vector<std::pair<int, int> > > unique_sets_of_pairs;
    std::vector<int> indices(all_pairs.size());
    std::iota(indices.begin(), indices.end(), 0); 
    // Generate all combinations of pairs and filter the valid ones
    if (balls.size()==4) {
        for (size_t i = 0; i < indices.size(); ++i) {
            for (size_t j = i + 1; j < indices.size(); ++j) {
                std::vector<std::pair<int, int>> candidate_set = {all_pairs[i], all_pairs[j]};
                if (isValidSet(candidate_set, balls)) {
                    std::sort(candidate_set.begin(), candidate_set.end());
                    unique_sets_of_pairs.insert(candidate_set);
                }   
            }   
        }   
    } else if (balls.size()==6) {
        for (size_t i = 0; i < indices.size(); ++i) {
            for (size_t j = i + 1; j < indices.size(); ++j) {
                for (size_t k = j + 1; k < indices.size(); ++k) {
                    std::vector<std::pair<int, int>> candidate_set = {all_pairs[i], all_pairs[j], all_pairs[k]};
                    if (isValidSet(candidate_set, balls)) {
                        std::sort(candidate_set.begin(), candidate_set.end());
                        unique_sets_of_pairs.insert(candidate_set);
                    }   
                }
            }
        }
    } else if (balls.size()==8) {
        for (size_t i = 0; i < indices.size(); ++i) {
            for (size_t j = i + 1; j < indices.size(); ++j) {
                for (size_t k = j + 1; k < indices.size(); ++k) {
                    for (size_t p = k + 1; p < indices.size(); ++p) {
                        std::vector<std::pair<int, int>> candidate_set = {all_pairs[i], all_pairs[j], all_pairs[k], all_pairs[p]};
                        if (isValidSet(candidate_set, balls)) {

                            std::sort(candidate_set.begin(), candidate_set.end());
                            unique_sets_of_pairs.insert(candidate_set);
                        }
                    }
                }
            }
        }
    } else if (balls.size()==10) {
        for (size_t i = 0; i < indices.size(); ++i) {
            for (size_t j = i + 1; j < indices.size(); ++j) {
                for (size_t k = j + 1; k < indices.size(); ++k) {
                    for (size_t p = k + 1; p < indices.size(); ++p) {
                        for (size_t q = p + 1; q < indices.size(); ++q) {
                            std::vector<std::pair<int, int>> candidate_set = {all_pairs[i], all_pairs[j], all_pairs[k], all_pairs[p], all_pairs[q]};
                            if (isValidSet(candidate_set, balls)) {
                                std::sort(candidate_set.begin(), candidate_set.end());
                                unique_sets_of_pairs.insert(candidate_set);
                            }
                        }
                    }
                }
            }
        }
    }
    int ii = int(uniqueSetsvec.at(num).size());
    int jj = 0;
    int npatterns = 0;
    // Step 4: Print the unique sets of pairs
    uniqueSetsvec.at(num).resize(uniqueSetsvec.at(num).size()+unique_sets_of_pairs.size());
    for (const auto& set_of_pairs : unique_sets_of_pairs) {
        uniqueSetsvec.at(num).at(ii).resize(set_of_pairs.size());
        jj = 0;
        for (const auto& p : set_of_pairs) {
            uniqueSetsvec.at(num).at(ii).at(jj).first = p.first;
            uniqueSetsvec.at(num).at(ii).at(jj).second = p.second;
            jj++;
        }
        npatterns++;
        ii++;
    }
    return;
}
// Function to check if a set of pairs is valid
bool WpMuonClassifyRecTool::isValidSet(const std::vector<std::pair<int, int> >& pairs, const std::vector<int>& balls) {
    std::map<int, int> counter;
    for (int i = 0; i < int(balls.size()); i++) {
        counter[balls.at(i)]++;
    }   
    std::map<int, int> count;
    for (const auto& p : pairs) {
        count[p.first]++;
        count[p.second]++;
    }   
    for (int i = 0; i < int(balls.size()); i++) {
        if (count[balls.at(i)] != counter[balls.at(i)]) return false;
    } 
    return true;
}