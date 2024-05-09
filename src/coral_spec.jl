@kwdef mutable struct CoralSpec
    f_groups = [
        "Tabular Acropora",
        "Corymbose Acropora",
        "Corymbose non-Acropora",
        "Small massives and encrusting",
        "Large massives"
    ]
    k_max::Float64 = 48671.09380827452
    settler_fracs::Vector{Float64} = [0.0384, 0.0315, 0.0681, 0.0363, 0.0363]
    initial_cover_fracs::Matrix{Float64} = [
        2.5113238193582336e-6 0.0002568268958279662 0.002730314784432863 0.011708198887993284 0.019926737221979646 0.013038774783840433 0.0004759662908532412;
        1.7178947155691395e-6 0.00017568485742776 0.0018676975481465923 0.008009103741809017 0.013631072308676492 0.008476250427529917 0.000768635514910798;
        7.059511430602026e-7 7.219588301626877e-5 0.0007675110744886387 0.0032912587076354303 0.005601548797045098 0.0033503122127778425 0.0004487784645370408;
        1.3805793685553956e-6 0.00014118844847369217 0.0015009678289966816 0.006436484894185198 0.010954557942400363 0.006557758903807724 0.0008718575221051794;
        1.526477814076502e-6 0.00015610912281304804 0.0016595888239322718 0.0071166871064378195 0.012112226245700793 0.007563408893376427 0.0006513630346739669
    ]
    linear_extension::Matrix{Float64} = [
        0.006094558 0.010718383 0.025514863 0.050798784 0.094509136 0.168505241 0.0;  # Tabular Acropora
        0.007685561 0.012208521 0.01864468 0.028229656 0.035293827 0.030042179 0.0;              # Corymbose Acropora
        0.001904555 0.003437468 0.006154666 0.009747701 0.017007948 0.029172889 0.0;      # Corymbose non-Acropora
        0.003180337 0.004738498 0.006837293 0.007105867 0.005810847 0.005810847 0.0;            # Small massives and encrusting
        0.001224784 0.00217702 0.003820976 0.00718781 0.012417243 0.020854626 0.0             # Large massives
    ]
    survival_rate::Matrix{Float64} = [
        0.859017851 0.858528906 0.857044217 0.856477498 0.856104353 0.855852241 0.855852241;    # Tabular Acropora
        0.865006527 0.87915437 0.892044073 0.905304164 0.915373252 0.925707536 0.925707536;     # Corymbose Acropora
        0.953069031 0.959152694 0.964460394 0.968306361 0.972598906 0.97621179 0.97621179;     # Corymbose non-Acropora
        0.869976692 0.938029324 0.977889252 0.987199004 0.99207702 0.996931548 0.996931548;     # Small massives and encrusting
        0.9782479 0.979496637 0.980850254 0.982178103 0.983568572 0.984667677 0.984667677       # Large massives
        # Large massives
    ]
    bins::Matrix{Float64} = [
        0.0 0.01 0.02 0.06 0.15 0.36 0.89 0.9;
        0.0 0.01 0.02 0.04 0.09 0.18 0.38 0.39;
        0.0 0.01 0.02 0.04 0.07 0.14 0.27 0.28;
        0.0 0.01 0.02 0.05 0.08 0.12 0.26 0.27;
        0.0 0.01 0.02 0.04 0.09 0.19 0.40 0.41
    ]
    # ADRIA
    # f_groups = [
    #     "Arborescent Acropora",
    #     "Tabular Acropora",
    #     "Corymbose Acropora",
    #     "Corymbose non-Acropora",
    #     "Small massives and encrusting",
    #     "Large massives"
    # ]
    # initial_cover_fracs::Matrix{Float64} = [
    #     1.1420056802592752e-6 0.00011679010552842463 0.0012415901799215276 0.0053242156716813045 0.009061534367300377 0.006145727188012592;
    #     2.5113238193582336e-6 0.0002568268958279662 0.002730314784432863 0.011708198887993284 0.019926737221979646 0.013514741074693675;
    #     1.7178947155691395e-6 0.00017568485742776 0.0018676975481465923 0.008009103741809017 0.013631072308676492 0.009244885942440715;
    #     7.059511430602026e-7 7.219588301626877e-5 0.0007675110744886387 0.0032912587076354303 0.005601548797045098 0.003799090677314883;
    #     1.3805793685553956e-6 0.00014118844847369217 0.0015009678289966816 0.006436484894185198 0.010954557942400363 0.007429616425912904;
    #     1.526477814076502e-6 0.00015610912281304804 0.0016595888239322718 0.0071166871064378195 0.012112226245700793 0.008214771928050394
    # ]
    # linear_extension::Matrix{Float64} = [
    #     0.01 0.03 0.03 0.044 0.044 0.044;       # Arborescent Acropora
    #     0.01 0.03 0.03 0.044 0.044 0.044;       # Tabular Acropora
    #     0.01 0.03 0.03 0.03 0.03 0.03;          # Corymbose Acropora
    #     0.01 0.024 0.024 0.024 0.024 0.024;     # Corymbose non-Acropora
    #     0.01 0.01 0.01 0.01 0.008 0.008;        # Small massives and encrusting
    #     0.01 0.01 0.01 0.01 0.012 0.012         # Large massives
    # ]
    # survival_rate::Matrix{Float64} = 1 .- [
    #     0.2 0.3 0.004 0.004 0.002 0.002;     # Arborescent Acropora
    #     0.2 0.3 0.190 0.190 0.098 0.098;     # Tabular Acropora
    #     0.2 0.3 0.172 0.172 0.088 0.088;     # Corymbose Acropora
    #     0.2 0.3 0.226 0.226 0.116 0.116;     # Corymbose non-Acropora
    #     0.2 0.3 0.040 0.026 0.020 0.020;     # Small massives and encrusting
    #     0.2 0.3 0.040 0.026 0.020 0.020     # Large massives
    # ]
end