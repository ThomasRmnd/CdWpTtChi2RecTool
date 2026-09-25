PKG(CdWpTtChi2RecTool

MODULE

FILES

    src/chi2/fht/FhtChi2.cpp
    src/chi2/tt/TtChi2.cpp

    src/cost/fht/FhtCostFunction.cpp
    src/cost/fht_tt/FhtTtCostFunction.cpp
    src/cost/tt/TtCostFunction.cpp

    src/estimator/Estimator.cpp
    src/estimator/ParamsModifierEstimator.cpp
    src/estimator/Pipeline.cpp
    src/estimator/fht/CorrectionMapLoopEstimator.cpp
    src/estimator/fht/FhtMinimizerEstimator.cpp
    src/estimator/fht/map/CorrectionMap.cpp
    src/estimator/fht/map/file/CorrectionFile.cpp
    src/estimator/fht_tt/FhtTtCorrMapEstimator.cpp
    src/estimator/fht/map/CorrectionMap.cpp
    src/estimator/fht/map/file/CorrectionFile.cpp
    src/estimator/optimizer/Optimizer.cpp
    src/estimator/optimizer/RootOptimizer.cpp
    src/estimator/tt/TtMinimizerEstimator.cpp
    src/estimator/tt/combinator/MaskHeightCartesianProdCombinator.cpp
    src/estimator/tt/converter/FuzeNeighborConverter.cpp

    # src/initializer/fht/ClusterBundleInitializer.cpp WIP
    src/initializer/fht/ClusterMaxChargeInitializer.cpp
    src/initializer/fht/WaterPhaseInitializer.cpp
    src/initializer/tt/MinMaxHeightInitializer.cpp

    src/strategy/Factory.cpp
    src/strategy/Registry.cpp
    src/strategy/Strategy.cpp

    src/transformer/CalibTimeTransformer.cpp
    src/transformer/EarlyLateFhtTransformer.cpp
    src/transformer/FhtChargeTholdTransformer.cpp
    src/transformer/PmtTypeTransformer.cpp
    src/transformer/Transformer.cpp
    src/transformer/TtCrossTalkTransformer.cpp
    src/transformer/WaterPhaseTransformer.cpp

    src/trigger/TriggerTimeCorrelation.cpp

    # src/utils/DensityPeaksClustering.cpp WIP
    src/utils/Geometry.cpp

DEPENDS

    GeometryLib
    Minuit2
    RecTools

)