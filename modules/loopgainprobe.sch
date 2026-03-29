v {xschem version=3.4.8RC file_version=1.3}
G {}
K {}
V {}
S {}
F {}
E {}
T {* The loopgainprobe subcircuit needs to be named X999 (numeric node name chosen for compatibility with KiCad)

* The probe subcircuit is inserted as follows, the a_node/b_node (input/output) order is unimportant
* X999 a_node b_node loopgainprobe

* The following .control section performs the loop analysis and plots the results
* .control 
* run 
* alter i.X999.Ii acmag=1 
* alter v.X999.Vi acmag=0 
* run
* gnuplot magfile db(tian_loop()) 
* gnuplot phfile 180*cph(tian_loop())/pi 
* .endc

* The loopgainprobe can be disabled by setting both i.X999.Ii and v.X999.Vi to zero

* tian_loop.lib was translated from Frank Wiedmann's LTspice LoopGainProbe.asc into the ngspice dialect
* For more information see: https://sites.google.com/site/frankwiedmann/loopgain

* Michael Tian, V. Visvanathan, Jeffrey Hantgan, and Kenneth Kundert,
*    "Striving for Small-Signal Stability", IEEE Circuits and Devices Magazine,
*     vol. 17, no. 1, pp. 31-41, January 2001.

.subckt loopgainprobe a b
Ii 0 x DC 0 AC 0
Vi x a DC 0 AC 1
Vnodebuffer b x 0
.ends loopgainprobe

.func tian_loop() \{1/(1-1/(2*(ac1.I(v.X999.Vi)*ac2.V(X999.x)-ac1.V(X999.x)*ac2.I(v.X999.Vi))+ac1.V(X999.x)+ac2.I(v.X999.Vi)))\}
} 240 -510 0 0 0.4 0.4 {}
T {.control
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
  meas ac gm_val find mag when phase=0

  let ugbw = ugbw_freq
  let pm = pm_val
  let gm = gm_val

  write stb.raw mag phase ugbw pm gm LoopGain
.endc} -580 0 0 0 0.4 0.4 {}
N -50 -80 50 -80 {lab=x}
N 0 -80 -0 -30 {lab=x}
N -150 -80 -110 -80 {lab=A}
N 110 -80 150 -80 {lab=B}
N 0 30 0 80 {lab=gnd}
C {isource.sym} 0 0 2 1 {name=Ii value="DC 0 AC 0"}
C {iopin.sym} 0 80 1 0 {name=p1 lab=gnd}
C {vsource.sym} -80 -80 1 1 {name=Vi value="DC 0 AC 0" savecurrent=false}
C {vsource.sym} 80 -80 1 1 {name=Vnodebuffer value=0 savecurrent=false}
C {iopin.sym} -150 -80 0 1 {name=p2 lab=A}
C {iopin.sym} 150 -80 0 0 {name=p3 lab=B}
C {lab_wire.sym} -10 -80 0 0 {name=p4 sig_type=std_logic lab=x}
