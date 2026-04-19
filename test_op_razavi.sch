v {xschem version=3.4.8RC file_version=1.3}
G {}
K {}
V {}
S {}
F {}
E {}
B 2 600 -540 1400 -140 {flags=graph
y1=0.053
y2=34
ypos1=0
ypos2=2
divy=5
subdivy=1
unity=1
x1=0
x2=12
divx=5
subdivx=8
xlabmag=1.0
ylabmag=1.0
dataset=-1
unitx=1
logx=1
logy=0
sim_type=ac
rawfile=$netlist_dir/test_op_razavi.raw
color=4
node=out}
N -80 -240 80 -240 {lab=cmfb1}
N -120 -210 -120 -30 {lab=n1}
N 120 -210 120 -30 {lab=p1}
N -120 70 120 70 {lab=s}
N -120 30 -120 70 {lab=s}
N 120 30 120 70 {lab=s}
N -0 70 -0 210 {lab=s}
N -320 180 320 180 {lab=cmfb2}
N -360 -300 360 -300 {lab=vdd}
N 360 -300 360 -270 {lab=vdd}
N 360 -270 360 -240 {lab=vdd}
N 120 -300 120 -240 {lab=vdd}
N -120 -300 -120 -240 {lab=vdd}
N -360 -300 -360 -240 {lab=vdd}
N -0 240 -0 300 {lab=0}
N -360 -210 -360 150 {lab=voutp}
N -360 180 -360 300 {lab=0}
N 360 -210 360 150 {lab=voutn}
N -360 300 360 300 {lab=0}
N 360 180 360 300 {lab=0}
N -600 270 -600 300 {lab=0}
N -180 -240 -180 -190 {lab=n1}
N 180 -240 180 -190 {lab=p1}
N -30 -180 30 -180 {lab=cmfb1}
N -120 -180 -90 -180 {lab=n1}
N 90 -180 120 -180 {lab=p1}
N 0 -240 0 -180 {lab=cmfb1}
N -360 120 -330 120 {lab=voutp}
N -270 120 -230 120 {lab=cmfb2}
N -230 120 -230 180 {lab=cmfb2}
N 240 120 270 120 {lab=cmfb2}
N 240 120 240 180 {lab=cmfb2}
N 330 120 360 120 {lab=voutn}
N -600 150 -600 210 {lab=#net1}
N -600 180 -540 180 {lab=#net1}
N -540 180 -540 240 {lab=#net1}
N -720 -300 -720 90 {lab=vdd}
N -720 150 -720 300 {lab=0}
N -720 300 -600 300 {lab=0}
N -120 0 -90 0 {lab=s}
N -90 0 -90 70 {lab=s}
N 90 -0 120 -0 {lab=s}
N 90 -0 90 70 {lab=s}
N 120 -190 180 -190 {lab=p1}
N 180 -240 320 -240 {lab=p1}
N -360 -190 -340 -190 {lab=voutp}
N -280 -190 -270 -190 {lab=#net2}
N 340 -190 360 -190 {lab=voutn}
N 270 -190 280 -190 {lab=#net3}
N 190 -190 210 -190 {lab=p1}
N -180 -190 -120 -190 {lab=n1}
N -320 -240 -180 -240 {lab=n1}
N -210 -190 -180 -190 {lab=n1}
N 180 -190 190 -190 {lab=p1}
N -560 240 -40 240 {lab=#net1}
N -600 300 -360 300 {lab=0}
N -720 -300 -360 -300 {lab=vdd}
N -240 -0 -160 -0 {lab=vip}
N 160 0 240 0 {lab=vin}
N -600 -300 -600 90 {lab=vdd}
N -960 270 -960 300 {lab=0}
N -960 150 -960 210 {lab=VCM}
N -600 240 -600 270 {lab=0}
N -1080 270 -1080 300 {lab=0}
N -1080 100 -1080 210 {lab=in}
N -1080 100 -1000 100 {lab=in}
N -1080 -20 -1080 100 {lab=in}
N -1080 -20 -880 -20 {lab=in}
N -960 180 -840 180 {lab=VCM}
N -840 30 -840 180 {lab=VCM}
N -840 -120 -840 -30 {lab=vip}
N -960 20 -960 90 {lab=vin}
N 600 -20 680 -20 {lab=voutp}
N 600 20 680 20 {lab=voutn}
N 720 -80 720 -30 {lab=out}
N 720 -80 790 -80 {lab=out}
N 790 -80 790 -30 {lab=out}
N 790 30 790 60 {lab=0}
N 720 60 790 60 {lab=0}
N 720 30 720 120 {lab=0}
C {sky130_fd_pr/nfet_01v8_lvt.sym} -140 0 0 0 {name=M1
L=0.15
W=48
nf=12 mult=1
model=nfet_01v8_lvt
spiceprefix=X
}
C {sky130_fd_pr/nfet_01v8_lvt.sym} 140 0 0 1 {name=M2
L=0.15
W=48
nf=12 mult=1
model=nfet_01v8_lvt
spiceprefix=X
}
C {sky130_fd_pr/nfet_01v8.sym} -340 180 0 1 {name=M7
L=0.2
W=24  
nf=6 mult=1
model=nfet_01v8
spiceprefix=X
}
C {sky130_fd_pr/nfet_01v8.sym} 340 180 0 0 {name=M8
L=0.2
W=24  
nf=6 mult=1
model=nfet_01v8
spiceprefix=X
}
C {sky130_fd_pr/pfet_01v8.sym} 100 -240 0 0 {name=M4
L=0.2
W=48
nf=12 mult=1
model=pfet_01v8
spiceprefix=X
}
C {sky130_fd_pr/pfet_01v8.sym} -100 -240 0 1 {name=M3
L=0.2
W=48
nf=12 mult=1
model=pfet_01v8
spiceprefix=X
}
C {sky130_fd_pr/pfet_01v8.sym} -340 -240 0 1 {name=M5
L=0.2
W=96
nf=24 mult=1
model=pfet_01v8
spiceprefix=X
}
C {sky130_fd_pr/pfet_01v8.sym} 340 -240 0 0 {name=M6
L=0.2
W=96
nf=24 mult=1
model=pfet_01v8
spiceprefix=X
}
C {sky130_fd_pr/res_high_po_0p35.sym} -60 -180 1 0 {name=R6
L=40
model=res_high_po_0p35
spiceprefix=X
 mult=1}
C {sky130_fd_pr/res_high_po_0p35.sym} 60 -180 1 0 {name=R1
L=40
model=res_high_po_0p35
spiceprefix=X
 mult=1}
C {sky130_fd_pr/res_high_po_0p35.sym} -300 120 1 1 {name=R2
L=40
model=res_high_po_0p35
spiceprefix=X
 mult=1}
C {sky130_fd_pr/res_high_po_0p35.sym} 300 120 1 1 {name=R3
L=40
model=res_high_po_0p35
spiceprefix=X
 mult=1}
C {isource.sym} -600 120 0 0 {name=I0 value=10u}
C {vsource.sym} -720 120 0 0 {name=V1 value=1.8 savecurrent=true}
C {capa.sym} -240 -190 1 0 {name=C1
m=1
value=200f
footprint=1206
device="ceramic capacitor"}
C {res.sym} -310 -190 1 0 {name=R4
value=400
footprint=1206
device=resistor
m=1}
C {capa.sym} 240 -190 3 1 {name=C2
m=1
value=200f
footprint=1206
device="ceramic capacitor"}
C {res.sym} 310 -190 3 1 {name=R5
value=400
footprint=1206
device=resistor
m=1}
C {gnd.sym} -720 300 0 0 {name=l1 lab=0}
C {lab_wire.sym} -690 -300 0 0 {name=p1 sig_type=std_logic lab=vdd}
C {lab_wire.sym} -120 -120 0 0 {name=p2 sig_type=std_logic lab=n1}
C {lab_wire.sym} 120 -120 0 0 {name=p3 sig_type=std_logic lab=p1}
C {lab_wire.sym} -360 -120 0 0 {name=p4 sig_type=std_logic lab=voutp}
C {lab_wire.sym} 360 -120 0 0 {name=p5 sig_type=std_logic lab=voutn}
C {lab_wire.sym} -30 -240 0 0 {name=p6 sig_type=std_logic lab=cmfb1}
C {lab_wire.sym} -120 180 0 0 {name=p7 sig_type=std_logic lab=cmfb2}
C {lab_pin.sym} -240 0 0 0 {name=p8 sig_type=std_logic lab=vip}
C {lab_pin.sym} 240 0 0 1 {name=p10 sig_type=std_logic lab=vin}
C {lab_wire.sym} 0 70 0 0 {name=p9 sig_type=std_logic lab=s}
C {lab_wire.sym} -960 20 0 0 {name=p11 sig_type=std_logic lab=vin}
C {gnd.sym} -960 300 0 0 {name=l3 lab=0}
C {vsource.sym} -960 240 0 0 {name=V3 value=0.9 savecurrent=false}
C {lab_wire.sym} -840 -120 0 0 {name=p12 sig_type=std_logic lab=vip}
C {code.sym} -1120 -260 0 0 {
name=TT_MODELS
only_toplevel=true
format="tcleval( @value )"
value="
** opencircuitdesign pdks install
.lib $::SKYWATER_MODELS/sky130.lib.spice tt
"
spice_ignore=false
      }
C {devices/simulator_commands_shown.sym} -1410 -280 0 0 {name=COMMANDS1
simulator=ngspice
only_toplevel=false 
value="
.include test_op_razavi.save
.option reltol=1e-5
+  abstol=1e-14 savecurrents
.control
  save all
  op
  remzerovec 
  write test_op_razavi.raw
  set appendwrite
  ac dec 10 1 1e12
  remzerovec
  write test_op_razavi.raw
#  tran 0.1n 100n
#  write test_op_razavi.raw
.endc
"}
C {sky130_fd_pr/annotate_fet_params.sym} -220 30 0 0 {name=annot1 ref=M1}
C {devices/launcher.sym} -880 -240 0 0 {name=h1
descr="Annotate OP" 
tclcommand="set show_hidden_texts 1; xschem annotate_op"
}
C {sky130_fd_pr/nfet_01v8.sym} -20 240 0 0 {name=M11
L=0.5
W=16  
nf=4 mult=1
model=nfet_01v8
spiceprefix=X
}
C {sky130_fd_pr/nfet_01v8.sym} -580 240 0 1 {name=M12
L=0.5
W=4  
nf=1 mult=1
model=nfet_01v8
spiceprefix=X
}
C {gnd.sym} -300 140 0 0 {name=l4 lab=0}
C {gnd.sym} 300 140 0 0 {name=l5 lab=0}
C {gnd.sym} -60 -200 2 1 {name=l6 lab=0}
C {gnd.sym} 60 -200 2 1 {name=l7 lab=0}
C {sky130_fd_pr/annotate_fet_params.sym} 130 40 0 0 {name=annot2 ref=M2}
C {sky130_fd_pr/annotate_fet_params.sym} 420 -260 0 0 {name=annot3 ref=M6
}
C {sky130_fd_pr/annotate_fet_params.sym} 420 180 0 0 {name=annot4 ref=M8}
C {sky130_fd_pr/annotate_fet_params.sym} 170 -150 0 0 {name=annot5 ref=M4
}
C {sky130_fd_pr/annotate_fet_params.sym} 60 200 0 0 {name=annot6 ref=M11}
C {vcvs.sym} -960 120 0 0 {name=E1 value=-0.5}
C {vcvs.sym} -840 0 0 0 {name=E2 value=0.5}
C {gnd.sym} -1080 300 0 0 {name=l2 lab=0}
C {vsource.sym} -1080 240 0 0 {name=V2 value="dc 0 ac 1" savecurrent=false}
C {gnd.sym} -1000 140 1 0 {name=l8 lab=0}
C {gnd.sym} -880 20 1 0 {name=l9 lab=0}
C {lab_wire.sym} -1080 -20 0 0 {name=p13 sig_type=std_logic lab=in}
C {lab_wire.sym} -960 200 0 0 {name=p14 sig_type=std_logic lab=VCM}
C {vcvs.sym} 720 0 0 0 {name=E3 value=-0.5}
C {gnd.sym} 720 120 0 0 {name=l10 lab=0}
C {lab_wire.sym} 610 -20 0 0 {name=p15 sig_type=std_logic lab=voutp}
C {lab_wire.sym} 610 20 0 0 {name=p16 sig_type=std_logic lab=voutn}
C {lab_wire.sym} 720 -80 0 0 {name=p17 sig_type=std_logic lab=out}
C {res.sym} 790 0 0 0 {name=R7
value=1M
footprint=1206
device=resistor
m=1}
