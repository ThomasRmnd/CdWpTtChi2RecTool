PKG(CdWpTtChi2RecTool

    src/estimator/Estimator.cpp
    src/estimator/ParamsModifierEstimator.cpp
    src/estimator/Pipeline.cpp
    src/estimator/fht/FhtMinimizerEstimator.cpp
    src/estimator/fht/map/file/CorrectionFile.cpp
    src/estimator/optimizer/NelderMeadOptimizer.cpp
    src/estimator/optimizer/Optimizer.cpp
    src/estimator/optimizer/RootOptimizer.cpp
    src/estimator/tt/TtMinimizerEstimator.cpp
    src/estimator/tt/combinator/MaskHeightCartesianProdCombinator.cpp
    src/estimator/tt/converter/FuzeNeighborConverter.cpp

    src/initializer/fht/ClusterBundleInitializer.cpp
    src/initializer/fht/ClusterMaxChargeInitializer.cpp
    src/initializer/fht/WaterPhaseInitializer.cpp
    src/initializer/tt/MinMaxHeightInitializer.cpp

    src/strategy/Factory.cpp
    src/strategy/Registry.cpp
    src/strategy/Strategy.cpp

    src/transformer/CalibrationTransformer.cpp
    src/transformer/EarlyLateFhtTransformer.cpp
    src/transformer/FhtChargeTholdTransformer.cpp
    src/transformer/PmtLocTypeTransformer.cpp
    src/transformer/Transformer.cpp
    src/transformer/TtCrossTalkTransformer.cpp
    src/transformer/WaterPhaseTransformer.cpp
    src/transformer/WpGeomTimeTransformer.cpp

    utils/DensityPeaksClustering.cpp
    utils/Geometry.cpp

DEPENDS
        
    Minuit2
    RecTools

)