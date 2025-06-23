PKG(CdWpTtChi2RecTool

    src/estimator/Estimator.cpp

    src/initializer/fht/ClusterBundleInitializer.cpp
    src/initializer/fht/ClusterMaxChargeInitializer.cpp
    src/initializer/fht/WaterPhaseInitializer.cpp
    src/initializer/tt/MinMaxHeightInitializer.cpp

    src/transformer/CalibrationTransformer.cpp
    src/transformer/EarlyLateFhtTransformer.cpp
    src/transformer/FhtChargeTholdTransformer.cpp
    src/transformer/PmtLocTypeTransformer.cpp
    src/transformer/Transformer.cpp

DEPENDS
        
    Minuit2
    RecTools

)