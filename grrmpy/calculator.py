def pfp_calculator(model_version="v3.0.0",calc_mode="crystal_plus_d3"):
    """PFP calculator

    Parameters:
    
    model_version: str
        modelのバージョン
    calc_mode: str
        | **'crystal'**: 
        |    結晶系(Hubbard補正あり)
        | **'crystal_plus_d3'**: 
        |   結晶系(Hubbard補正あり)+DFT-D3補正
        | **'crystal_u0'**:
        |   結晶系(Hubbard補正なし)
        | **'crystal_u0_plus_d3'**: 
        |   結晶系(Hubbard補正なし)+DFT-D3補正
        | **'molecule'**:
        |   分子系
    """
    from pfp_api_client.pfp.calculators.ase_calculator import ASECalculator
    from pfp_api_client.pfp.estimator import Estimator, EstimatorCalcMode
    if calc_mode == "crystal":
        calc_mode = EstimatorCalcMode.CRYSTAL
    elif calc_mode == "crystal_plus_d3":
        calc_mode = EstimatorCalcMode.CRYSTAL_PLUS_D3
    elif calc_mode == "crystal_u0":
        calc_mode = EstimatorCalcMode.CRYSTAL_U0
    elif calc_mode == "crystal_u0_plus_d3":
        calc_mode = EstimatorCalcMode.CRYSTAL_U0_PLUS_D3
    elif calc_mode == "molecule":
        calc_mode = EstimatorCalcMode.MOLECULE
    else:
        calc_mode = calc_mode 
    estimator = Estimator(calc_mode=calc_mode,model_version=model_version)
    calc = ASECalculator(estimator)
    return calc


