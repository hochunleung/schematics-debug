v {xschem version=3.4.8RC file_version=1.3}
G {}
K {}
V {}
S {}
F {}
E {}
B 2 170 -400 970 0 {flags=graph
y1=0
y2=2
ypos1=0
ypos2=2
divy=5
subdivy=1
unity=1
x1=0
x2=1.8
divx=5
subdivx=1
xlabmag=1.0
ylabmag=1.0
node="in
out"
color="7 4"
dataset=-1
unitx=1
logx=0
logy=0
sim_type=dc
rawfile=$netlist_dir/untitled.raw
autoload=1}
B 2 170 0 970 400 {flags=graph
y1=0
y2=2
ypos1=0
ypos2=2
divy=5
subdivy=1
unity=1
x1=0
x2=5e-06
divx=5
subdivx=1
xlabmag=1.0
ylabmag=1.0
node="in
out"
color="8 12"
dataset=-1
unitx=1
logx=0
logy=0
sim_type=tran
rawfile=$netlist_dir/untitled.raw
autoload=1}
N 0 -50 0 50 {lab=out}
N -0 110 -0 140 {lab=GND}
N -0 80 20 80 {lab=GND}
N 20 80 20 140 {lab=GND}
N 0 140 20 140 {lab=GND}
N -180 140 0 140 {lab=GND}
N -180 -140 -0 -140 {lab=vcc}
N -0 -140 -0 -110 {lab=vcc}
N -220 -140 -220 40 {lab=vcc}
N -220 100 -220 140 {lab=GND}
N -140 100 -140 140 {lab=GND}
N -60 80 -40 80 {lab=in}
N -60 -80 -60 80 {lab=in}
N -60 -80 -40 -80 {lab=in}
N -140 0 -140 40 {lab=in}
N -100 -0 -60 -0 {lab=in}
N -140 0 -100 0 {lab=in}
N -220 -140 -180 -140 {lab=vcc}
N -220 140 -180 140 {lab=GND}
C {sky130_fd_pr/nfet_01v8.sym} -20 80 0 0 {name=M2
L=0.15
W=1  
nf=1 mult=1
model=nfet_01v8
spiceprefix=X
}
C {sky130_fd_pr/pfet_01v8.sym} -20 -80 0 0 {name=M1
L=0.15
W=1
nf=1 mult=1
model=pfet_01v8
spiceprefix=X
}
C {gnd.sym} 0 140 0 0 {name=l1 lab=GND}
C {vsource.sym} -140 70 0 0 {name=VIN value="DC 0.9 sin(0.8 0.1 1meg 0 0 0)" savecurrent=true}
C {vsource.sym} -220 70 0 0 {name=VCC value=1.8 savecurrent=true}
C {lab_pin.sym} -140 0 0 0 {name=p1 sig_type=std_logic lab=in}
C {lab_pin.sym} -150 -140 0 0 {name=p2 sig_type=std_logic lab=vcc}
C {lab_pin.sym} 0 0 0 0 {name=p3 sig_type=std_logic lab=out}
C {code.sym} -450 50 0 0 {
name=TT_MODELS
only_toplevel=true
format="tcleval( @value )"
value="
** opencircuitdesign pdks install
.lib $::SKYWATER_MODELS/sky130.lib.spice tt
"
spice_ignore=false
      }
C {sky130_fd_pr/annotate_fet_params.sym} 60 50 0 0 {name=annot1 ref=M2}
C {sky130_fd_pr/annotate_fet_params.sym} 60 -110 0 0 {name=annot2 ref=M1}
C {sky130_fd_pr/annotate_fet_params.sym} 60 -110 0 0 {name=annot3 ref=M1}
C {devices/simulator_commands_shown.sym} -530 -310 0 0 {name=COMMANDS1
simulator=ngspice
only_toplevel=false 
value="
.include untitled.save
.option reltol=1e-5
+  abstol=1e-14 savecurrents
.control
  save all
  op
  remzerovec 
  write untitled.raw
  set appendwrite
  tran 0.1u 5u
  remzerovec
  write untitled.raw
  set appendwrite
  dc VIN 0 1.8 0.1
  write untitled.raw
.endc
"}
C {devices/launcher.sym} -160 -280 0 0 {name=h3
descr="Netlist & sim" 
tclcommand="xschem netlist; xschem simulate"}
C {devices/launcher.sym} -160 -190 0 0 {name=h2 
descr="Annotate OP
Load/unload TRAN waveforms" 
tclcommand="
xschem raw_read $netlist_dir/untitled.raw tran;
xschem raw_read $netlist_dir/untitled.raw dc;
set show_hidden_texts 1; xschem annotate_op
"
}
C {devices/launcher.sym} -160 -240 0 0 {name=h1
descr="Annotate OP" 
tclcommand="set show_hidden_texts 1; xschem annotate_op"
}
C {devices/launcher.sym} -160 -320 0 0 {name=h5
descr="Generate .save lines" 
tclcommand="write_data [sky130_save_fet_params] $netlist_dir/[file rootname [file tail [xschem get current_name]]].save
textwindow $netlist_dir/[file rootname [file tail [xschem get current_name]]].save
"
}
