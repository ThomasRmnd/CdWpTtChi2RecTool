PKG(CdWpTtChi2RecTool

    src/estimator/Estimator.cpp
    src/estimator/ParamsModifierEstimator.cpp
    src/estimator/Pipeline.cpp

    src/initializer/fht/ClusterBundleInitializer.cpp
    src/initializer/fht/ClusterMaxChargeInitializer.cpp
    src/initializer/fht/WaterPhaseInitializer.cpp
    src/initializer/tt/MinMaxHeightInitializer.cpp

    src/transformer/CalibrationTransformer.cpp
    src/transformer/EarlyLateFhtTransformer.cpp
    src/transformer/FhtChargeTholdTransformer.cpp
    src/transformer/PmtLocTypeTransformer.cpp
    src/transformer/Transformer.cpp
    src/transformer/TtCrossTalkTransformer.cpp
    src/transformer/WaterPhaseTransformer.cpp
    src/transformer/WpGeomtryTimeTransformer.cpp

DEPENDS
        
    Minuit2
    RecTools

)