v {xschem version=3.4.8RC file_version=1.3}
G {}
K {}
V {}
S {}
F {}
E {}
N -0 -140 -0 -110 {lab=out}
N 80 -140 80 -110 {lab=out}
N 160 -140 160 -110 {lab=out}
N 0 -140 160 -140 {lab=out}
N 0 -0 80 -0 {lab=vcm}
N 80 -0 160 -0 {lab=vcm}
N -120 -100 -40 -100 {lab=vip}
N -120 -60 -40 -60 {lab=vin}
N 160 -140 200 -140 {lab=out}
N -120 0 -0 0 {lab=vcm}
N 0 -50 0 -0 {lab=vcm}
N 80 -50 80 -0 {lab=vcm}
N 160 -50 160 -0 {lab=vcm}
C {vccs.sym} 0 -80 0 0 {name=G1 value=\{gm\}}
C {res.sym} 80 -80 0 0 {name=R1
value=\{ro1\}
footprint=1206
device=resistor
m=1}
C {capa.sym} 160 -80 0 0 {name=C1
m=1
value=\{c1\}
footprint=1206
device="ceramic capacitor"}
C {ipin.sym} -120 -100 0 0 {name=p1 lab=vip}
C {ipin.sym} -120 -60 0 0 {name=p2 lab=vin}
C {opin.sym} 200 -140 0 0 {name=p4 lab=out}
C {ipin.sym} -120 0 0 0 {name=p5 lab=vcm}
