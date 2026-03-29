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
N -0 -140 -0 -110 {lab=#net1}
N 80 -140 80 -110 {lab=#net1}
N 160 -140 160 -110 {lab=#net1}
N 0 -140 160 -140 {lab=#net1}
N 0 -0 80 -0 {lab=vcm}
N 80 -0 160 -0 {lab=vcm}
N 0 110 0 140 {lab=#net2}
N 80 110 80 140 {lab=#net2}
N 160 110 160 140 {lab=#net2}
N -0 140 160 140 {lab=#net2}
N -120 -100 -40 -100 {lab=vip}
N -80 -100 -80 60 {lab=vip}
N -80 60 -40 60 {lab=vip}
N -120 -60 -40 -60 {lab=vin}
N -60 -60 -60 100 {lab=vin}
N -60 100 -40 100 {lab=vin}
N -120 0 -0 0 {lab=vcm}
N 360 -50 360 50 {lab=vcm}
N 440 -50 440 50 {lab=vcm}
N 360 -140 360 -110 {lab=outp}
N 440 -140 440 -110 {lab=outp}
N 360 0 440 0 {lab=vcm}
N 360 110 360 140 {lab=outn}
N 440 110 440 140 {lab=outn}
N 160 -140 200 -140 {lab=#net1}
N 200 -140 200 -100 {lab=#net1}
N 200 -100 240 -100 {lab=#net1}
N 160 -0 360 -0 {lab=vcm}
N 160 140 200 140 {lab=#net2}
N 200 100 200 140 {lab=#net2}
N 200 100 240 100 {lab=#net2}
N 240 0 240 60 {lab=vcm}
N 240 -60 240 -0 {lab=vcm}
N 280 -50 280 -0 {lab=vcm}
N 280 0 280 50 {lab=vcm}
N 280 -140 280 -110 {lab=outp}
N 280 -140 510 -140 {lab=outp}
N 280 110 280 140 {lab=outn}
N 280 140 510 140 {lab=outn}
N 310 -220 360 -220 {lab=outp}
N 360 -220 360 -140 {lab=outp}
N 190 -220 250 -220 {lab=#net3}
N 80 -220 130 -220 {lab=#net1}
N 80 -220 80 -140 {lab=#net1}
N 310 220 360 220 {lab=outn}
N 190 220 250 220 {lab=#net4}
N 80 220 130 220 {lab=#net2}
N 360 140 360 220 {lab=outn}
N 80 140 80 220 {lab=#net2}
C {vccs.sym} 0 -80 0 0 {name=G1 value=\{gm1/2\}}
C {vccs.sym} 0 80 0 0 {name=G2 value=\{gm1/2\}}
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
value=\{co1\}
footprint=1206
device="ceramic capacitor"}
C {capa.sym} 160 80 0 0 {name=C2
m=1
value=\{co1\}
footprint=1206
device="ceramic capacitor"}
C {ipin.sym} -120 -100 0 0 {name=p1 lab=vip}
C {ipin.sym} -120 -60 0 0 {name=p2 lab=vin}
C {opin.sym} 510 -140 0 0 {name=p3 lab=outp}
C {opin.sym} 510 140 0 0 {name=p4 lab=outn}
C {ipin.sym} -120 0 0 0 {name=p5 lab=vcm}
C {vccs.sym} 280 -80 0 0 {name=G3 value=\{gm2\}}
C {vccs.sym} 280 80 0 0 {name=G4 value=\{gm2\}}
C {res.sym} 360 -80 0 0 {name=R3
value=\{ro2\}
footprint=1206
device=resistor
m=1}
C {res.sym} 360 80 0 0 {name=R4
value=\{ro2\}
footprint=1206
device=resistor
m=1}
C {capa.sym} 440 -80 0 0 {name=C3
m=1
value=\{co2\}
footprint=1206
device="ceramic capacitor"}
C {capa.sym} 440 80 0 0 {name=C4
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
C {res.sym} 160 220 3 1 {name=R6
value=\{rz\}
footprint=1206
device=resistor
m=1}
C {capa.sym} 280 220 3 1 {name=C6
m=1
value=\{cc\}
footprint=1206
device="ceramic capacitor"}
