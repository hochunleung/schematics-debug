v {xschem version=3.4.8RC file_version=1.3}
G {}
K {}
V {}
S {}
F {}
E {}
N -0 -140 -0 -110 {lab=#net1}
N 80 -140 80 -110 {lab=#net1}
N 160 -140 160 -110 {lab=#net1}
N 0 -140 160 -140 {lab=#net1}
N 0 -0 80 -0 {lab=vcm}
N 80 -0 160 -0 {lab=vcm}
N -120 -100 -40 -100 {lab=vip}
N -120 -60 -40 -60 {lab=vin}
N -120 0 -0 0 {lab=vcm}
N 360 -140 360 -110 {lab=out}
N 440 -140 440 -110 {lab=out}
N 360 0 440 0 {lab=vcm}
N 160 -140 200 -140 {lab=#net1}
N 200 -140 200 -100 {lab=#net1}
N 200 -100 240 -100 {lab=#net1}
N 160 -0 360 -0 {lab=vcm}
N 240 -60 240 -0 {lab=vcm}
N 280 -50 280 -0 {lab=vcm}
N 280 -140 280 -110 {lab=out}
N 280 -140 510 -140 {lab=out}
N 310 -220 360 -220 {lab=out}
N 360 -220 360 -140 {lab=out}
N 190 -220 250 -220 {lab=#net2}
N 80 -220 130 -220 {lab=#net1}
N 80 -220 80 -140 {lab=#net1}
N 0 -50 -0 -0 {lab=vcm}
N 80 -50 80 -0 {lab=vcm}
N 160 -50 160 -0 {lab=vcm}
N 360 -50 360 -0 {lab=vcm}
N 440 -50 440 0 {lab=vcm}
C {vccs.sym} 0 -80 0 0 {name=G1 value=\{gm1\}}
C {res.sym} 80 -80 0 0 {name=R1
value=\{ro1\}
footprint=1206
device=resistor
m=1}
C {capa.sym} 160 -80 0 0 {name=C1
m=1
value=\{co1\}
footprint=1206
device="ceramic capacitor"}
C {ipin.sym} -120 -100 0 0 {name=p1 lab=vip}
C {ipin.sym} -120 -60 0 0 {name=p2 lab=vin}
C {opin.sym} 510 -140 0 0 {name=p3 lab=out}
C {ipin.sym} -120 0 0 0 {name=p5 lab=vcm}
C {vccs.sym} 280 -80 0 0 {name=G3 value=\{gm2\}}
C {res.sym} 360 -80 0 0 {name=R3
value=\{ro2\}
footprint=1206
device=resistor
m=1}
C {capa.sym} 440 -80 0 0 {name=C3
m=1
value=\{co2\}
footprint=1206
device="ceramic capacitor"}
C {res.sym} 160 -220 3 0 {name=R5
value=\{rz\}
footprint=1206
device=resistor
m=1}
C {capa.sym} 280 -220 3 0 {name=C5
m=1
value=\{cc\}
footprint=1206
device="ceramic capacitor"}
