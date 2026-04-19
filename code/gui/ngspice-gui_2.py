import customtkinter as ctk
import sys
import subprocess
import threading
import os

# 設定全局主題與風格
ctk.set_appearance_mode("dark")  # 深色模式更能體現散光效果
ctk.set_default_color_theme("blue") # 也可以自定義主題文件

class NgspiceGUI(ctk.CTk):
    def __init__(self, netlist):
        super().__init__()

        self.title(f"ngspice-gui Pro - {netlist}")
        self.geometry("1100x800") # 稍微加大視窗以容納更多設置

        # 設定字體：使用你習慣的微軟正黑體 Light，展現優雅感
        self.main_font = ("Microsoft JhengHei Light", 14)
        self.header_font = ("Microsoft JhengHei", 20, "bold")

        # 處理檔名與路徑
        self.netlist_name = os.path.basename(netlist)
        self.base_name = os.path.splitext(self.netlist_name)[0]
        
        self.abs_netlist_path = os.path.abspath(netlist)
        self.work_dir = os.path.dirname(self.abs_netlist_path)
        self.env = os.environ.copy()
        self.env["NGSPICE_MEAS_RAWFILE"] = "1" # 減少內部緩衝

        # 模組化配置設計：未來可輕鬆移至 JSON 文件
        self.analysis_config = [
            {
                "id": "op",
                "title": "OP Analysis",
                "fields": [],
                "default_file": f"{self.base_name}.raw",
                "cmd_template": "op"
            },
            {
                "id": "dc",
                "title": "DC Analysis",
                "fields": [("Source", "v1"), ("Start", "0"), ("Stop", "1.8"), ("Step", "0.01")],
                "default_file": f"{self.base_name}_dc.raw",
                "cmd_template": "dc {Source} {Start} {Stop} {Step}"
            },
            {
                "id": "ac",
                "title": "AC Analysis",
                "fields": [("Type", "dec"), ("Points", "10"), ("fStart", "1"), ("fStop", "100meg")],
                "default_file": f"{self.base_name}_ac.raw",
                "cmd_template": "ac {Type} {Points} {fStart} {fStop}"
            },
            {
                "id": "tran",
                "title": "Transient Analysis",
                "fields": [("Tstep", "1n"), ("Tstop", "1u")],
                "default_file": f"{self.base_name}_tran.raw",
                "cmd_template": "tran {Tstep} {Tstop}"
            }
        ]

        # 儲存各個分析區域的變數，用於預覽同步
        self.analysis_vars = {}

        # 啟動 ngspice 背景進程
        self.process = subprocess.Popen(
            # 使用 stdbuf -oL 強制 stdout 為行緩衝，確保橫幅訊息不會被緩衝
            ['stdbuf', '-oL', '-eL', 'ngspice', '-i', self.abs_netlist_path],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, # 重要：將錯誤與正常輸出合併，解決順序問題
            text=True,
            bufsize=1,                # 行緩衝
            shell=False,
            cwd=self.work_dir,
            env=self.env
        )

        self.setup_ui()

        # 啟動異步讀取
        threading.Thread(target=self.read_output, daemon=True).start()

    def setup_ui(self):
        # 配置網格權重
        self.grid_columnconfigure(0, weight=0) # 左側設置區固定寬度
        self.grid_columnconfigure(1, weight=1) # 右側 Log 區自動拉伸
        self.grid_rowconfigure(0, weight=1)

        # 創建左側滾動設置區
        self.settings_scroll = ctk.CTkScrollableFrame(self, width=450, corner_radius=0, fg_color="transparent")
        self.settings_scroll.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)

        # 動態生成仿真介面
        self.analysis_data = {}
        for config in self.analysis_config:
            self.analysis_data[config['id']] = self.create_analysis_group(
                config['title'],
                config['fields'],
                config['default_file'],
                config['cmd_template']
            )

        # 5. 輸出日誌區域 (移到右側)
        self.log_area = ctk.CTkTextbox(
            self, 
            corner_radius=15, 
            border_width=2, 
            border_color="#3B3B3B",
            font=("Consolas", 12),
            fg_color="#0D0D0D"
        )
        self.log_area.grid(row=0, column=1, padx=20, pady=20, sticky="nsew")

    def create_analysis_group(self, title, fields, default_file, cmd_template):
        """輔助函式：創建帶有標題、輸入框、檔名和按鈕的組合框"""
        frame = ctk.CTkFrame(self.settings_scroll, corner_radius=20, border_width=2, border_color="#3B3B3B")
        frame.pack(fill="x", padx=10, pady=10)

        label = ctk.CTkLabel(frame, text=title, font=self.header_font)
        label.pack(pady=(10, 5))

        # 參數網格
        input_grid = ctk.CTkFrame(frame, fg_color="transparent")
        input_grid.pack(padx=10, pady=5)

        field_vars = {}
        entry_kwargs = {
            "width": 80,
            "corner_radius": 8,
            "border_width": 2,
            "font": ("Consolas", 12)
        }
        
        # 更新預覽的函式
        def update_preview(*args):
            # 1. 取得核心仿真命令
            params = {k: v.get() for k, v in field_vars.items()}
            sim_cmd = cmd_template.format(**params)
            
            all_cmds = [sim_cmd] # 仿真命令放在最前面
            if save_var.get():
                # 2. 配置命令 (按照要求在仿真之後，remzerovec 不使用 set)
                if rem_var.get(): all_cmds.append("remzerovec")
                if app_var.get(): all_cmds.append("set appendwrite")
                else: all_cmds.append("unset appendwrite")
                
                # 3. 寫入命令
                all_cmds.append(f"write {file_var.get()}")

            preview_entry.configure(state="normal")
            preview_entry.delete(0, "end")
            preview_entry.insert(0, " | ".join(all_cmds))
            preview_entry.configure(state="readonly")

        # 繪製參數欄位
        for i, (f_label, default) in enumerate(fields):
            lbl = ctk.CTkLabel(input_grid, text=f_label, font=self.main_font)
            lbl.grid(row=i // 2, column=(i % 2) * 2, padx=5, pady=5, sticky="e")
            
            var = ctk.StringVar(value=default)
            entry = ctk.CTkEntry(input_grid, textvariable=var, **entry_kwargs)
            entry.grid(row=i // 2, column=(i % 2) * 2 + 1, padx=5, pady=5)
            field_vars[f_label] = var
            var.trace_add("write", update_preview)

        # 存檔設置區域 (聯動邏輯)
        save_var = ctk.BooleanVar(value=True)
        save_check = ctk.CTkCheckBox(frame, text="Save result?", variable=save_var, font=self.main_font)
        save_check.pack(pady=5)

        file_settings_frame = ctk.CTkFrame(frame, fg_color="transparent")
        file_settings_frame.pack(fill="x", padx=10, pady=5)

        # 檔名輸入
        file_entry_frame = ctk.CTkFrame(file_settings_frame, fg_color="transparent")
        file_entry_frame.pack(fill="x", pady=2)
        ctk.CTkLabel(file_entry_frame, text="Raw File:", font=self.main_font).pack(side="left")
        
        file_var = ctk.StringVar(value=default_file)
        file_entry = ctk.CTkEntry(file_entry_frame, textvariable=file_var, corner_radius=8, border_width=2, font=("Consolas", 12))
        file_entry.pack(side="left", padx=5, fill="x", expand=True)
        file_var.trace_add("write", update_preview)

        # remzerovec 與 append 勾選框
        opts_frame = ctk.CTkFrame(file_settings_frame, fg_color="transparent")
        opts_frame.pack(fill="x", pady=2)
        rem_var = ctk.BooleanVar(value=True)
        rem_check = ctk.CTkCheckBox(opts_frame, text="remzerovec", variable=rem_var, font=("Consolas", 11), width=100)
        rem_check.pack(side="left", padx=5)
        rem_var.trace_add("write", update_preview)

        app_var = ctk.BooleanVar(value=False)
        app_check = ctk.CTkCheckBox(opts_frame, text="append", variable=app_var, font=("Consolas", 11), width=100)
        app_check.pack(side="left", padx=5)
        app_var.trace_add("write", update_preview)

        # 控制顯示/隱藏的函式
        def toggle_file_settings(*args):
            if save_var.get():
                file_settings_frame.pack(fill="x", padx=10, pady=5, before=run_btn)
            else:
                file_settings_frame.pack_forget()
            update_preview()
        
        save_var.trace_add("write", toggle_file_settings)

        # 命令預覽 (放在按鈕上方)
        preview_frame = ctk.CTkFrame(frame, fg_color="transparent")
        preview_frame.pack(fill="x", padx=20, pady=5)
        preview_entry = ctk.CTkEntry(preview_frame, font=("Consolas", 11), state="readonly", border_width=0, fg_color="#1A1A1A", text_color="#AAAAAA")
        preview_entry.pack(fill="x", expand=True)

        update_preview() # 初始化預覽

        # 執行按鈕
        run_btn = ctk.CTkButton(
            frame, 
            text=f"RUN {title.split()[0]}", 
            command=lambda: self.execute_simulation(title, fields, field_vars, save_var, file_var, rem_var, app_var, cmd_template),
            font=self.main_font,
            corner_radius=10,
            fg_color="#1f538d",
            hover_color="#2c73c3"
        )
        run_btn.pack(pady=10)

        return {}

    def send_command(self, cmd):
        if self.process.poll() is None:
            try:
                # 如果指令包含 write，ngspice 會覆蓋文件
                self.process.stdin.write(cmd + "\n")
                self.process.stdin.flush()
                # 過濾掉 set/unset 命令，讓 Log 保持乾淨，只顯示核心指令
                if not any(x in cmd for x in ["set ", "unset ", "write ", "remzerovec"]):
                    self.write_log(f"\n> {cmd}\n")
            except Exception:
                pass

    def execute_simulation(self, title, fields, field_vars, save_var, file_var, rem_var, app_var, cmd_template):
        """通用的模擬執行邏輯，包含 write 命令處理"""
        # 1. 取得並執行仿真命令 (確保它在最前面)
        params = {k: v.get() for k, v in field_vars.items()}
        sim_cmd = cmd_template.format(**params)
        self.send_command(sim_cmd)

        if save_var.get():
            # 2. 發送配置命令 (放在仿真之後)
            if rem_var.get(): 
                self.send_command("remzerovec") # 直接發送，不使用 set
            
            if app_var.get(): self.send_command("set appendwrite")
            else: self.send_command("unset appendwrite")

            # 3. 發送寫入命令
            self.send_command(f"write {file_var.get()}")

    def read_output(self):
        # 使用 iter 確保即時讀取每一行
        while True:
            line = self.process.stdout.readline()
            if line:
                # 過濾掉 \r (Carriage Return)，解決方塊顯示問題
                clean_line = line.replace('\r', '')
                self.write_log(clean_line)
            elif self.process.poll() is not None:
                # 進程結束且無輸出時退出
                break

    def write_log(self, text):
        # 使用 after 確保在主線程更新 UI
        self.after(0, self._update_log_area, text)

    def _update_log_area(self, text):
        self.log_area.insert("end", text)
        self.log_area.see("end")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python ngspice-gui.py <netlist>")
        sys.exit(1)
    
    app = NgspiceGUI(sys.argv[1])
    app.mainloop()