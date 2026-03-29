v {xschem version=3.4.8RC file_version=1.3}
G {}
K {}
V {}
S {}
F {}
E {}
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
C {vcvs.sym} -180 50 0 0 {name=E1 value=0.5}
C {vcvs.sym} -120 110 0 0 {name=E2 value=-0.5}
C {vsource.sym} -120 260 0 0 {name=V1 value=0.9 savecurrent=false}
C {vsource.sym} -260 130 0 0 {name=V2 value="0 ac 1" savecurrent=false}
C {gnd.sym} -120 320 0 0 {name=l1 lab=GND}
C {lab_pin.sym} -120 200 0 0 {name=p1 sig_type=std_logic lab=vcm}
C {lab_pin.sym} -70 30 0 0 {name=p2 sig_type=std_logic lab=vcm}
C {lab_pin.sym} -180 -10 0 0 {name=p4 sig_type=std_logic lab=vip}
C {lab_pin.sym} -120 10 0 0 {name=p5 sig_type=std_logic lab=vin}
C {devices/simulator_commands_shown.sym} -560 60 0 0 {name=COMMANDS1
simulator=ngspice
only_toplevel=false 
value="
.include test.save
.option reltol=1e-5
+  abstol=1e-14 savecurrents
.control
  save all
  op
  remzerovec 
  write test.raw
.endc
"}
C {gnd.sym} -260 180 0 0 {name=l2 lab=GND}
