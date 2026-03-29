v {xschem version=3.4.8RC file_version=1.3}
G {}
K {}
V {}
S {}
F {}
E {}
N 0 -50 -0 50 {lab=vcm}
N 80 -50 80 50 {lab=vcm}
N 160 -50 160 50 {lab=vcm}
N -0 -140 -0 -110 {lab=outn}
N 80 -140 80 -110 {lab=outn}
N 160 -140 160 -110 {lab=outn}
N 0 -140 160 -140 {lab=outn}
N 0 -0 80 -0 {lab=vcm}
N 80 -0 160 -0 {lab=vcm}
N 0 110 0 140 {lab=outp}
N 80 110 80 140 {lab=outp}
N 160 110 160 140 {lab=outp}
N -0 140 160 140 {lab=outp}
N -120 -100 -40 -100 {lab=vip}
N -80 -100 -80 60 {lab=vip}
N -80 60 -40 60 {lab=vip}
N -120 -60 -40 -60 {lab=vin}
N -60 -60 -60 100 {lab=vin}
N -60 100 -40 100 {lab=vin}
N 160 -140 200 -140 {lab=outn}
N 160 140 200 140 {lab=outp}
N -120 0 -0 0 {lab=vcm}
C {vccs.sym} 0 -80 0 0 {name=G1 value=\{gm/2\}}
C {vccs.sym} 0 80 0 0 {name=G2 value=\{gm/2\}}
C {res.sym} 80 -80 0 0 {name=R1
value=\{ro1\}
footprint=1206
device=resistor
m=1}
C {res.sym} 80 80 0 0 {name=R2
value=\{ro1\}
footprint=1206
device=resistor
m=1}
C {capa.sym} 160 -80 0 0 {name=C1
m=1
value=\{c1\}
footprint=1206
device="ceramic capacitor"}
C {capa.sym} 160 80 0 0 {name=C2
m=1
value=\{c1\}
footprint=1206
device="ceramic capacitor"}
C {ipin.sym} -120 -100 0 0 {name=p1 lab=vip}
C {ipin.sym} -120 -60 0 0 {name=p2 lab=vin}
C {opin.sym} 200 140 0 0 {name=p3 lab=outp}
C {opin.sym} 200 -140 0 0 {name=p4 lab=outn}
C {ipin.sym} -120 0 0 0 {name=p5 lab=vcm}
