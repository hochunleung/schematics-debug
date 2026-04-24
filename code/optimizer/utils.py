import json

def load_circuit_config(config_path='config.json'):
    """
    從 JSON 配置文件中載入並解析電路仿真參數。
    """
    with open(config_path, 'r', encoding='utf-8') as f:
        config = json.load(f)

    # 提取技術參數
    tech_params = config.get('tech_params', {})
    vdd_val = tech_params.get('vdd', 0.9)
    vcm_val = tech_params.get('vcm', 0.45)
    stb_probe = tech_params.get('stb_probe', False)
    is_include_stb = str(stb_probe).strip().lower() in ('true', '1', 't', 'y', 'yes')

    # 建立映射字典
    voltage_map = {
        'vdd': vdd_val,
        'vcm': vcm_val,
        'gnd': 0.0
    }

    # 動態生成 fixed_voltages
    fixed_node_types = config.get('fixed_voltages', {})
    fixed_voltages = {}
    for node, val in fixed_node_types.items():
        if isinstance(val, (int, float)):
            fixed_voltages[node] = float(val)
        elif isinstance(val, str):
            val_lower = val.lower()
            if val_lower in voltage_map:
                fixed_voltages[node] = voltage_map[val_lower]
            else:
                try:
                    fixed_voltages[node] = float(val)
                except ValueError:
                    raise ValueError(f"⚠️ 配置錯誤: 無法解析節點 '{node}' 的電壓設定 '{val}'")
        else:
            raise TypeError(f"⚠️ 配置錯誤: 節點 '{node}' 的電壓類型不支援 '{type(val)}'")
    
    saved_currents = config.get('saved_currents', [])
    solver_settings = config.get('solver_settings', {})
    #ignored_isources = solver_settings.get('ignored_isources', [])
    #kcl_tricks = solver_settings.get('kcl_tricks', {})

    return {
        'vdd_val': vdd_val,
        'vcm_val': vcm_val,
        'is_include_stb': is_include_stb,
        'fixed_voltages': fixed_voltages,
        'saved_currents': saved_currents,
        'solver_settings': solver_settings,
        #'ignored_isources': ignored_isources,
        #'kcl_tricks': kcl_tricks,
        'raw_config': config
    }