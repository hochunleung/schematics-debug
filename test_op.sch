v {xschem version=3.4.8RC file_version=1.3}
G {}
K {}
V {}
S {}
F {}
E {}
B 2 260 -300 1060 100 {flags=graph
y1=-33
y2=60
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
B 2 1060 -300 1860 100 {flags=graph,unlocked
y1=0.024
y2=1000
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
node=vout
color=4
dataset=-1
unitx=1
logx=1
logy=0
autoload=1
sim_type=ac}
B 2 260 100 1060 500 {flags=graph
y1=-180
y2=-1.8e-05
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
sim_type=ac}
N -120 10 -120 80 {lab=vin}
N -120 10 -60 10 {lab=vin}
N -180 -10 -180 20 {lab=vip}
N -180 -10 -60 -10 {lab=vip}
N -180 80 -180 180 {lab=vcm}
N -180 180 -120 180 {lab=vcm}
N -120 140 -120 180 {lab=vcm}
N -260 30 -260 100 {lab=#net1}
N -260 30 -220 30 {lab=#net1}
N -260 90 -150 90 {lab=#net1}
N -220 70 -220 180 {lab=GND}
N -260 180 -220 180 {lab=GND}
N -260 160 -260 180 {lab=GND}
N -220 130 -160 130 {lab=GND}
N -120 180 -120 230 {lab=vcm}
N -120 290 -120 320 {lab=GND}
N -90 30 -60 30 {lab=vcm}
N 60 -10 100 -10 {lab=voutp}
N 60 10 80 10 {lab=voutn}
N 160 10 160 30 {lab=voutn}
N 160 30 180 30 {lab=voutn}
N 220 40 220 60 {lab=GND}
N 220 -40 220 -20 {lab=vout}
N 80 10 160 10 {lab=voutn}
N 100 -10 180 -10 {lab=voutp}
C {vcvs.sym} -180 50 0 0 {name=E1 value=0.5}
C {vcvs.sym} -120 110 0 0 {name=E2 value=-0.5}
C {vsource.sym} -120 260 0 0 {name=V1 value=0.9 savecurrent=false}
C {vsource.sym} -260 130 0 0 {name=V2 value="0 ac 1" savecurrent=false}
C {gnd.sym} -120 320 0 0 {name=l1 lab=GND}
C {lab_pin.sym} -120 200 0 0 {name=p1 sig_type=std_logic lab=vcm}
C {lab_pin.sym} -70 30 0 0 {name=p2 sig_type=std_logic lab=vcm}
C {vcvs.sym} 220 10 0 0 {name=E3 value=1}
C {gnd.sym} 220 60 0 0 {name=l2 lab=GND}
C {lab_pin.sym} 220 -40 0 0 {name=p3 sig_type=std_logic lab=vout}
C {lab_pin.sym} -180 -10 0 0 {name=p4 sig_type=std_logic lab=vip}
C {lab_pin.sym} -120 10 0 0 {name=p5 sig_type=std_logic lab=vin}
C {lab_pin.sym} 120 -10 0 0 {name=p6 sig_type=std_logic lab=voutp}
C {lab_pin.sym} 120 10 0 0 {name=p7 sig_type=std_logic lab=voutn}
C {devices/simulator_commands_shown.sym} -560 60 0 0 {name=COMMANDS1
simulator=ngspice
only_toplevel=false 
value="
.include test_op.save
.option reltol=1e-5
+  abstol=1e-14 savecurrents
.control
  save all
  op
  remzerovec 
  write test_op.raw
  set appendwrite
  ac dec 10 1 1e10
  remzerovec 
  write test_op.raw
.endc
"}
C {gnd.sym} -260 180 0 0 {name=l3 lab=GND}
C {devices/launcher.sym} -10 140 0 0 {name=h2 
descr="Annotate OP
Load/unload AC waveforms" 
tclcommand="
xschem raw_read $netlist_dir/untitled.raw ac;
set show_hidden_texts 1; xschem annotate_op
"
}
C {/home/tester/proj/modules/ideal_op_diff_2stage.sym} 0 0 0 0 {name=x1 gm1=2.5m ro1=8k co1=100f gm2=10m ro2=5k co2=2p cc=100f rz=200}
