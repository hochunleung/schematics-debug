v {xschem version=3.4.8RC file_version=1.3}
G {}
K {}
V {}
S {}
F {}
E {}
B 2 300 -280 1100 120 {flags=graph
y1=-24.06
y2=15.64
ypos1=0
ypos2=2
divy=5
subdivy=1
unity=1
x1=0
x2=10
divx=5
subdivx=8
xlabmag=1.0
ylabmag=1.0
node="vout db20()"
color=4
dataset=-1
unitx=1
logx=1
logy=0
autoload=1
sim_type=ac
}
B 2 300 120 1100 520 {flags=graph
y1=-180
y2=-1.8e-08
ypos1=0
ypos2=2
divy=5
subdivy=1
unity=1
x1=0
x2=10
divx=5
subdivx=8
xlabmag=1.0
ylabmag=1.0
node=ph(vout)
color=7
dataset=-1
unitx=1
logx=1
logy=0
sim_type=ac
}
B 2 1100 -280 1900 120 {flags=graph
ypos1=0
ypos2=2
divy=5
subdivy=1
unity=1
x1=0
x2=10
divx=5
subdivx=8
xlabmag=1.0
ylabmag=1.0
node="loopgain db20()"
color=4
dataset=-1
unitx=1
logx=1
logy=0
rawfile=$netlist_dir/stb.raw
sim_type=ac
y1=-20
y2=60}
B 2 1100 120 1900 520 {flags=graph
y1=0
y2=180
ypos1=0
ypos2=2
divy=5
subdivy=1
unity=1
x1=0
x2=10
divx=5
subdivx=8
xlabmag=1.0
ylabmag=1.0
node=phase
color=7
dataset=-1
unitx=1
logx=1
logy=0
rawfile=$netlist_dir/stb.raw
sim_type=ac}
N -120 10 -60 10 {lab=vout}
N -180 -10 -150 -10 {lab=vip}
N -180 -10 -180 60 {lab=vip}
N -180 190 -180 220 {lab=GND}
N -90 30 -60 30 {lab=vcm}
N 60 0 100 0 {lab=#net1}
N -180 110 -180 130 {lab=vcm}
N -150 -10 -60 -10 {lab=vip}
N -120 -100 -120 10 {lab=vout}
N -120 -100 100 -100 {lab=vout}
N 220 -100 220 0 {lab=vout}
N 100 -100 220 -100 {lab=vout}
N 160 50 160 80 {lab=GND}
C {vsource.sym} -180 160 0 0 {name=V1 value=0.9 savecurrent=false}
C {vsource.sym} -180 80 0 0 {name=V2 value="0 ac 1" savecurrent=false}
C {gnd.sym} -180 220 0 0 {name=l1 lab=GND}
C {lab_pin.sym} -180 120 0 0 {name=p1 sig_type=std_logic lab=vcm}
C {lab_pin.sym} -90 30 0 0 {name=p2 sig_type=std_logic lab=vcm}
C {lab_pin.sym} -180 -10 0 0 {name=p4 sig_type=std_logic lab=vip}
C {lab_pin.sym} 220 -100 0 1 {name=p6 sig_type=std_logic lab=vout}
C {devices/launcher.sym} -10 140 0 0 {name=h2 
descr="Annotate OP
Load/unload AC waveforms" 
tclcommand="
xschem raw_read $netlist_dir/test_op_single.raw ac;
set show_hidden_texts 1; xschem annotate_op
"
}
C {/home/tester/proj/modules/ideal_op_single_2stage.sym} 0 0 0 0 {name=x1 gm1=2.5m ro1=8k co1=100f gm2=10m ro2=5k co2=2p cc=100f rz=200}
C {gnd.sym} 160 80 0 0 {name=l2 lab=GND}
C {devices/simulator_commands_shown.sym} -560 -120 0 0 {name=COMMANDS1
simulator=ngspice
only_toplevel=false 
value="
.include test_op_single_fb.save
.option reltol=1e-5
+  abstol=1e-14 savecurrents
.control
  save all
  op
  remzerovec 
  write test_op_single_fb.raw
  set appendwrite
  ac dec 10 1 1e10
  remzerovec 
  write test_op_single_fb.raw
  alter V2 acmag=0
.endc
"}
C {/home/tester/proj/modules/loopgainprobe.sym} 160 0 0 0 {name=X999}
C {simulator_commands_shown.sym} -560 300 0 0 {
name=COMMANDS_STB
simulator=ngspice
only_toplevel=false 
value="
.control
  shell rm -f stb.raw
  alter i.X999.Ii acmag=0
  alter v.X999.Vi acmag=1
  ac dec 10 1 1e10
  set v_run = $curplot

  alter i.X999.Ii acmag=1
  alter v.X999.Vi acmag=0
  ac dec 10 1 1e10
  set i_run = $curplot
  setplot $i_run

  let LoopGain = 1/(1-1/(2*(\{$v_run\}.I(v.X999.Vi)*\{$i_run\}.V(X999.x)-\{$v_run\}.V(X999.x)*\{$i_run\}.I(v.X999.Vi))+\{$v_run\}.V(X999.x)+\{$i_run\}.I(v.X999.Vi)))
  let mag = db(LoopGain)
  let phase = 180*cph(LoopGain)/pi

  meas ac ugbw_freq when mag=0
  meas ac pm_val find phase when mag=0

  let ugbw = ugbw_freq
  let pm = pm_val

  write stb.raw mag phase ugbw pm LoopGain
.endc
"
      }
