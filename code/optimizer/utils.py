import json

def load_circuit_config(config_path='config.json'):
    """
    從 JSON 配置文件中載入並解析電路仿真參數。
    """
    with open(config_path, 'r', encoding='utf-8') as f:
        config = json.load(f)

    # 提取技術參數
    dc = config.get('DC', {})
    if 'tech_params' in dc:
        tech_params = dc.get('tech_params', {})
        vdd_val = tech_params.get('vdd', 0.9)
        vcm_val = tech_params.get('vcm', 0.45)

        # 建立映射字典
        voltage_map = {
            'vdd': vdd_val,
            'vcm': vcm_val,
            'gnd': 0.0
        }

        # 動態生成 fixed_voltages
        fixed_node_types = dc.get('fixed_voltages', {})
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
        dc['fixed_voltages'] = fixed_voltages
        stb_probe = tech_params.get('stb_probe', False)
        is_include_stb = str(stb_probe).strip().lower() in ('true', '1', 't', 'y', 'yes')
        dc['tech_params']['is_include_stb'] = is_include_stb

    # saved_currents = config.get('saved_currents', [])
    # saved_op = config.get('saved_op', [])
    # solver_settings = config.get('solver_settings', {})
    #ignored_isources = solver_settings.get('ignored_isources', [])
    #kcl_tricks = solver_settings.get('kcl_tricks', {})

    stb = config.get('STB', {})
    is_run_CM = stb.get('is_run_CM', False)
    is_run_CM = str(is_run_CM).strip().lower() in ('true', '1', 't', 'y', 'yes')
    stb['is_run_CM'] = is_run_CM

    hd = config.get('Distor', {})
    noise = config.get('noise', {})

    return {
        'DC': dc,
        'STB': stb,
        'hd': hd,
        'noise': noise
    }