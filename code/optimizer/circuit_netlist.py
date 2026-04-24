import re
import numpy as np
from pprint import pprint

class CircuitNetlist:
    dict_model_to_type = {
        'resistor': 'R',
        'capacitor': 'C',
        'nch_mac': 'NMOS',
        'pch_mac': 'PMOS',
        'nch_lvt_mac': 'NMOS',
        'vsource': 'VSOURCE',
        'isource': 'ISOURCE',
        'pcccs': 'PCCCS',
        'vcvs': 'VCVS'
    }
    dict_param_interest = {
        'NMOS': ['w', 'l'],
        'PMOS': ['w', 'l'],
        'R': ['r'],
        'C': ['c']
    }
    def __init__(self, netlist_path=None, dialect='spectre'):
        self.elements = []  # 儲存所有元件物件
        self.dialect = dialect.lower() # 'spectre' 或 'ngspice'
        self.nodes = set()  # 所有節點名稱
        self.node_to_idx = {}
        self.config = {
            'fixed_voltages': {},
            'input_nodes': [],
            'output_nodes': [],
            'gnd_node': [] #'0'
        }
        if netlist_path:
            self.parse_file(netlist_path)

    def parse_file(self, filepath):
        """解析類 SPICE 語法的 netlist 檔案"""
        with open(filepath, 'r') as f:
            raw_lines = f.readlines()

        # 1. 處理分行描述 (Line Continuation)
        # 1. 處理分行描述 (Line Continuation)
        lines = []
        if self.dialect == 'ngspice':
            for line in raw_lines:
                s_line = line.strip()
                if not s_line: continue
                if s_line.startswith('+'):
                    if lines:
                        lines[-1] += " " + s_line[1:].strip()
                else:
                    lines.append(s_line)
        else:  # spectre 模式
            curr_line = ""
            for line in raw_lines:
                s_line = line.strip()
                if not s_line: continue
                if s_line.endswith('\\'):
                    curr_line += s_line[:-1].strip() + " "
                else:
                    curr_line += s_line
                    lines.append(curr_line)
                    curr_line = ""

        # 2. 忽略仿真命令 (Simulation Commands)
        in_control_block = False # Ngspice 專用
        for line in lines:
            # 根據 dialect 處理註釋
            comment_chars = ('*', '#', '//') if self.dialect == 'spectre' else ('*', '#')
            if any(line.startswith(c) for c in comment_chars):
                continue
            
            # Ngspice .control 塊處理
            if self.dialect == 'ngspice' and line.startswith('.'):
                if line.lower().startswith('.control'):
                    in_control_block = True
                elif line.lower().startswith('.endc'):
                    in_control_block = False
                # 忽略所有以 . 開頭的 Ngspice 指令 (如 .tran, .options, .temp, .lib, .save, .global)
                continue
            
            if self.dialect == 'ngspice' and in_control_block:
                # 忽略 .control 塊內的所有命令 (如 save, op, meas, plot, let, alter, ac, set, write, remzerovec)
                continue

            # if self.dialect == 'spectre':
            #     parts = line.split()
            #     if not parts: continue
            #     # Spectre 指令關鍵字
            #     spectre_cmds = {'simulator', 'global', 'parameters', 'include', 'library', 'saveoptions',
            #                     'dcOp', 'dcOpInfo', 'ac', 'tran', 'finalTimeOP', 'pss', 'modelParameter',
            #                     'element', 'outputParameter', 'designParamVals', 'primitives', 'subckts', 'options', 'save', 'statistics', 'noise', 'stb', 'model'}
            #     if parts[0].lower() in spectre_cmds:
            #         continue

            # 3. 提取元件名稱 (例如 VPRB, M1) 與初步過濾
            name_match = re.match(r'^(\S+)', line)
            if not name_match: continue
            
            # 4. 提取節點列表 (Spectre 核心篩選邏輯)
            nodes = []
            node_match = re.search(r'\((.*?)\)', line)
            
            if self.dialect == 'spectre' and not node_match:
                # 在 Spectre 語法中，器件實例必須有 (nodes)
                # 缺少括號的行通常是 parameters, options, save 或分析指令，直接跳過
                continue

            if node_match:
                nodes = node_match.group(1).strip().split()
                remainder = line[node_match.end():].strip()
            else:
                # 處理 Ngspice/SPICE 標準位置節點格式 (無括號)
                parts = line.split()
                type_char = parts[0][0].upper()
                # 在 SPICE 中，不同元件節點數量不同 (M=4, R/C/V/I=2, Q=3/4)
                node_count_map = {'M': 4, 'Q': 3, 'R': 2, 'C': 2, 'L': 2, 'V': 2, 'I': 2}
                n_count = node_count_map.get(type_char, 2) # 預設取 2
                
                if len(parts) > n_count:
                    nodes = parts[1:n_count+1]
                    remainder = " ".join(parts[n_count+1:])
                else:
                    nodes = parts[1:]
            
            name = name_match.group(1).upper()
            type_char = name[0]

            # 5. 解析剩餘部分：Model Name 和 Parameters
            # remainder 可能為 "vsource dc=0 type=dc" 或 "nch_mac l=100n ..."
            rem_parts = re.split(r'\s+', remainder)
            model = ""
            param_start_idx = 0
            if rem_parts and '=' not in rem_parts[0] and not rem_parts[0].startswith('['):
                model = rem_parts[0]
                param_start_idx = 1
            
            params = {}
            param_str = " ".join(rem_parts[param_start_idx:])
            
            # 提取列表型參數 (例如 probes=[ vprb Vi ])
            list_params = re.findall(r'(\w+)\s*=\s*\[\s*(.*?)\s*\]', param_str)
            for k, v in list_params:
                params[k.lower()] = v.split()
            
            # 提取鍵值對參數 (例如 dc=0, l=100n)
            clean_param_str = re.sub(r'\w+\s*=\s*\[.*?\]', '', param_str)
            kv_params = re.findall(r'(\w+)\s*=\s*([^\s=]+)', clean_param_str)
            for k, v in kv_params:
                params[k.lower()] = self._parse_value(v)

            params_selected = {}            
            type = self.dict_model_to_type[model] if model in self.dict_model_to_type else model.upper()
            if type in self.dict_param_interest:
                for k in self.dict_param_interest[type]:
                    params_selected[k] = params[k]
            else:
                params_selected = params

            element = {
                'name': name,
                'type': type,
                'model': model,
                'nodes': nodes,
                'params': params_selected
            }
            self.elements.append(element)
            for node in nodes:
                if node != self.config['gnd_node']:
                    self.nodes.add(node)

    def _parse_value(self, val_str):
        """處理 SPICE 數值單位 (u, n, f, k, M)"""
        # 根據 Dialect 定義單位映射
        if self.dialect == 'spectre':
            # Spectre: m=milli, M=Mega
            unit_map = {'u': 1e-6, 'n': 1e-9, 'f': 1e-15, 'p': 1e-12, 'K': 1e3, 'k': 1e3, 'm': 1e-3, 'M': 1e6}
        else:
            # Ngspice: m=milli, M=milli, Meg=Mega
            unit_map = {'u': 1e-6, 'n': 1e-9, 'f': 1e-15, 'p': 1e-12, 'K': 1e3, 'k': 1e3, 'm': 1e-3, 'M': 1e-3, 'meg': 1e6}

        s = val_str.lower()
        if s.endswith('meg'):
            res = float(s[:-3]) * 1e6
            return float(f"{res:.12g}")
        
        # 檢查原始字串的最後一個字元 (區分 Spectre 的 m/M)
        last_char = val_str[-1]
        lookup_char = last_char if self.dialect == 'spectre' else last_char.lower()

        if lookup_char in unit_map:
            try:
                res = float(val_str[:-1]) * unit_map[lookup_char]
                return float(f"{res:.12g}")
            except ValueError:
                return val_str
        try:
            return float(s)
        except ValueError:
            return val_str

    def setup_indices(self, unknown_nodes=None):
        """
        建立節點到矩陣索引的映射。
        unknown_nodes: 如果指定，則只針對這些節點編號 (用於 DC 求解器)
        """
        target_nodes = unknown_nodes if unknown_nodes else sorted(list(self.nodes))
        self.node_to_idx = {node: i for i, node in enumerate(target_nodes)}
        return self.node_to_idx

    def update_element_params(self, param_map):
        """
        根據優化器的輸入動態更新元件參數。
        param_map: {'M1': {'w': 96e-6, 'l': 100e-9}, 'R34': {'r': 4000}}
        """
        for elem in self.elements:
            if elem['name'] in param_map:
                elem['params'].update(param_map[elem['name']])

    def get_elements_by_type(self, type_char):
        """篩選元件，例如獲取所有 'M' (MOSFET)"""
        return [e for e in self.elements if e['type'] == type_char.upper()]

    def get_element_by_name(self, name):
        for e in self.elements:
            if e['name'] == name.upper():
                return e
        return None

    def apply_config(self, config):
        """應用外部配置，指定輸入輸出節點等"""
        self.config.update(config)

if __name__ == "__main__":
    # 測試腳本
    #parser = CircuitNetlist('buf_razavi.net', dialect='spectre')
    parser = CircuitNetlist('buf_razavi_full.mdl', dialect='spectre')
    #parser = CircuitNetlist('inv.net', dialect='spectre')
    #parser = CircuitNetlist('buf_razavi_simp.net', dialect='spectre')
    #parser = CircuitNetlist('test_amp.spice', dialect='ngspice')

    # 模擬優化器更新參數
    # parser.update_element_params({
    #     'M1': {'w': 80e-6, 'l': 100e-9},
    #     'R34': {'r': 5000}
    # })
    
    # 建立索引
    #idx_map = parser.setup_indices()
    #print(f"Node Indices: {idx_map}")
    #pprint(idx_map)
    
    # 獲取所有 MOSFET 並查看其節點與寬度
    #for type in ['NMOS', 'PMOS', 'C', 'R', 'VSOURCE', 'ISOURCE', 'VCCS', 'PCCCS']:
    elements = parser.elements
    for e in elements:
        print(f"name: {e['name']}, type: {e['type']}, model: {e['model']}, nodes: {e['nodes']}, params: {e['params']}")
    
    pprint(parser.nodes)
    #pprint(parser.setup_indices())
    # 獲取 R38 的電阻值
    #r38 = parser.get_element_by_name('R38')
    #print(f"R38 Resistance: {r38['params'].get('r')}")