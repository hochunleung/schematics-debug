v {xschem version=3.4.8RC file_version=1.3}
G {}
K {}
V {}
S {}
F {}
E {}
N -50 0 50 0 {lab=x}
N 0 0 0 50 {lab=x}
N -150 0 -110 0 {lab=in1}
N 110 0 150 0 {lab=out1}
N 0 110 0 160 {lab=gnd}
N -180 -300 0 -300 {lab=in2}
N 0 -300 0 -250 {lab=in2}
N 0 -190 0 -100 {lab=out2}
N 0 -100 160 -100 {lab=out2}
N -90 -140 -0 -140 {lab=out2}
N -180 -140 -150 -140 {lab=gnd}
N -180 -140 -180 -100 {lab=gnd}
N -80 -240 -40 -240 {lab=in1}
N -80 -200 -40 -200 {lab=out1}
C {isource.sym} 0 80 2 1 {name=Ii value="DC 0 AC 0"}
C {iopin.sym} 0 160 1 0 {name=p1 lab=gnd}
C {vsource.sym} -80 0 1 1 {name=Vi value="DC 0 AC 0" savecurrent=false}
C {vsource.sym} 80 0 3 1 {name=Vprb value=0 savecurrent=false}
C {iopin.sym} -150 0 0 1 {name=p2 lab=in1}
C {iopin.sym} 150 0 0 0 {name=p3 lab=out1}
C {lab_wire.sym} -10 0 0 0 {name=p4 sig_type=std_logic lab=x}
C {vcvs.sym} 0 -220 0 0 {name=E1 value=-1}
C {lab_wire.sym} -60 -240 0 0 {name=p5 sig_type=std_logic lab=in1}
C {lab_wire.sym} -60 -200 0 0 {name=p6 sig_type=std_logic lab=out1}
C {cccs.sym} -120 -140 3 0 {name=F1 vnam="POLY(2) Vprb Vi" value="0 -1 -1"}
C {lab_wire.sym} -180 -100 0 0 {name=p7 sig_type=std_logic lab=gnd}
C {iopin.sym} -180 -300 0 1 {name=p8 lab=in2}
C {iopin.sym} 160 -100 0 0 {name=p9 lab=out2}
