//genesis
// kkit Version 11 flat dumpfile
 
// Saved on Tue Apr 30 10:00:23 2024
 
include kkit {argv 1}
 
FASTDT = 0.001
SIMDT = 0.001
CONTROLDT = 0.1
PLOTDT = 0.1
MAXTIME = 100
TRANSIENT_TIME = 2
VARIABLE_DT_FLAG = 0
DEFAULT_VOL = 6.5692e-20
VERSION = 11.0
setfield /file/modpath value ~/scripts/modules
kparms
 
//genesis

initdump -version 3 -ignoreorphans 1
simobjdump doqcsinfo filename accessname accesstype transcriber developer \
  citation species tissue cellcompartment methodology sources \
  model_implementation model_validation x y z
simobjdump table input output alloced step_mode stepsize x y z
simobjdump xtree path script namemode sizescale
simobjdump xcoredraw xmin xmax ymin ymax
simobjdump xtext editable
simobjdump xgraph xmin xmax ymin ymax overlay
simobjdump xplot pixflags script fg ysquish do_slope wy
simobjdump group xtree_fg_req xtree_textfg_req plotfield expanded movealone \
  link savename file version md5sum mod_save_flag x y z
simobjdump geometry size dim shape outside xtree_fg_req xtree_textfg_req x y \
  z
simobjdump kpool DiffConst CoInit Co n nInit mwt nMin vol slave_enable \
  geomname xtree_fg_req xtree_textfg_req x y z
simobjdump kreac kf kb notes xtree_fg_req xtree_textfg_req x y z
simobjdump kenz CoComplexInit CoComplex nComplexInit nComplex vol k1 k2 k3 \
  keepconc usecomplex notes xtree_fg_req xtree_textfg_req link x y z
simobjdump stim level1 width1 delay1 level2 width2 delay2 baselevel trig_time \
  trig_mode notes xtree_fg_req xtree_textfg_req is_running x y z
simobjdump xtab input output alloced step_mode stepsize notes editfunc \
  xtree_fg_req xtree_textfg_req baselevel last_x last_y is_running x y z
simobjdump kchan perm gmax Vm is_active use_nernst notes xtree_fg_req \
  xtree_textfg_req x y z
simobjdump transport input output alloced step_mode stepsize dt delay clock \
  kf xtree_fg_req xtree_textfg_req x y z
simobjdump proto x y z
simobjdump text str
simundump geometry /kinetics/geometry 0 6.5933e-20 3 sphere "" white black 0 \
  0 0
simundump geometry /kinetics/geometry[1] 0 6.5692e-20 3 sphere "" white black \
  0 0 0
simundump text /kinetics/notes 0 ""
call /kinetics/notes LOAD \
""
simundump group /kinetics/glu 0 blue green x 0 0 "" defaultfile defaultfile.g \
  0 0 0 19 5 0
simundump text /kinetics/glu/notes 0 Compartment
call /kinetics/glu/notes LOAD \
"Compartment"
simundump kpool /kinetics/glu/Ca 0 0.0 0.080003 0.080003 3.1649 3.1649 0 0 \
  39.56 0 /kinetics/geometry 62 yellow 20 3 0
simundump text /kinetics/glu/Ca/notes 0 Compartment
call /kinetics/glu/Ca/notes LOAD \
"Compartment"
simundump kpool /kinetics/glu/Ca_ext 0 0.0 0.080003 0.080003 3.1649 3.1649 0 \
  0 39.56 4 /kinetics/geometry 52 yellow 24 3 0
simundump text /kinetics/glu/Ca_ext/notes 0 ""
call /kinetics/glu/Ca_ext/notes LOAD \
""
simundump kpool /kinetics/glu/RR_pool 0 0.0 0.3 0.3 11.868 11.868 0 0 39.56 0 \
  /kinetics/geometry 1 yellow 18 1 0
simundump text /kinetics/glu/RR_pool/notes 0 ""
call /kinetics/glu/RR_pool/notes LOAD \
""
simundump kpool /kinetics/glu/vesicle_pool 0 0.0 1.289 1.289 50.992 50.992 0 \
  0 39.56 4 /kinetics/geometry 27 yellow 24 -1 0
simundump text /kinetics/glu/vesicle_pool/notes 0 ""
call /kinetics/glu/vesicle_pool/notes LOAD \
""
simundump kpool /kinetics/glu/glu 0 0.0 0 0 0 0 0 0 39.56 0 \
  /kinetics/geometry 7 yellow 24 1 0
simundump text /kinetics/glu/glu/notes 0 ""
call /kinetics/glu/glu/notes LOAD \
""
simundump kpool /kinetics/glu/Docked 0 0.0 0 0 0 0 0 0 39.56 0 \
  /kinetics/geometry 54 blue 22 1 0
simundump text /kinetics/glu/Docked/notes 0 ""
call /kinetics/glu/Docked/notes LOAD \
""
simundump kpool /kinetics/glu/Ca_RR 0 0.0 0 0 0 0 0 0 39.56 0 \
  /kinetics/geometry 51 blue 20 1 0
simundump text /kinetics/glu/Ca_RR/notes 0 ""
call /kinetics/glu/Ca_RR/notes LOAD \
""
simundump kreac /kinetics/glu/remove_Ca 0 0.37917 0.12639 "" white yellow 22 \
  4 0
simundump text /kinetics/glu/remove_Ca/notes 0 ""
call /kinetics/glu/remove_Ca/notes LOAD \
""
simundump kreac /kinetics/glu/remove 0 50000 0 "" white yellow 22 0 0
simundump text /kinetics/glu/remove/notes 0 ""
call /kinetics/glu/remove/notes LOAD \
""
simundump kreac /kinetics/glu/replenish_vesicle 0 2 2 "" white yellow 22 -2 0
simundump text /kinetics/glu/replenish_vesicle/notes 0 ""
call /kinetics/glu/replenish_vesicle/notes LOAD \
""
simundump kreac /kinetics/glu/vesicle_release 0 0.00063898 0 "" white yellow \
  23 2 0
simundump text /kinetics/glu/vesicle_release/notes 0 \
  "High cooperativity, 4 or higher. Several refs: McDargh and O-Shaughnessy, BioRxiv 2021 Voleti, Jaczynska, Rizo, eLife 2020 Chen.... Scheller, Cell 1999"
call /kinetics/glu/vesicle_release/notes LOAD \
"High cooperativity, 4 or higher. Several refs: McDargh and O-Shaughnessy, BioRxiv 2021 Voleti, Jaczynska, Rizo, eLife 2020 Chen.... Scheller, Cell 1999"
simundump kreac /kinetics/glu/Ca_bind_RR 0 0.0063898 2 "" white blue 19 2 0
simundump text /kinetics/glu/Ca_bind_RR/notes 0 ""
call /kinetics/glu/Ca_bind_RR/notes LOAD \
""
simundump kreac /kinetics/glu/docking 0 1500 0 "" white blue 21 2 0
simundump text /kinetics/glu/docking/notes 0 ""
call /kinetics/glu/docking/notes LOAD \
""
simundump kreac /kinetics/glu/undocking 0 200 0 "" white blue 20 0 0
simundump text /kinetics/glu/undocking/notes 0 ""
call /kinetics/glu/undocking/notes LOAD \
""
simundump kpool /kinetics/glu/Buffer 0 0 8 8 315.32 315.32 0 0 39.415 0 \
  /kinetics/geometry[1] 47 blue 16 3 0
simundump text /kinetics/glu/Buffer/notes 0 ""
call /kinetics/glu/Buffer/notes LOAD \
""
simundump kpool /kinetics/glu/Ca4.Buffer 0 0 0 0 0 0 0 0 39.415 0 \
  /kinetics/geometry[1] 50 blue 16 1 0
simundump text /kinetics/glu/Ca4.Buffer/notes 0 ""
call /kinetics/glu/Ca4.Buffer/notes LOAD \
""
simundump kreac /kinetics/glu/Ca_bind_buffer 0 4.1434e-06 1000 "" white blue \
  17 2 0
simundump text /kinetics/glu/Ca_bind_buffer/notes 0 ""
call /kinetics/glu/Ca_bind_buffer/notes LOAD \
""
simundump group /kinetics/GABA 0 blue green x 0 0 "" defaultfile \
  defaultfile.g 0 0 0 6 5 0
simundump text /kinetics/GABA/notes 0 Compartment
call /kinetics/GABA/notes LOAD \
"Compartment"
simundump kpool /kinetics/GABA/Ca 0 0.0 0.080003 0.080003 3.1649 3.1649 0 0 \
  39.56 0 /kinetics/geometry 62 yellow 7 3 0
simundump text /kinetics/GABA/Ca/notes 0 ""
call /kinetics/GABA/Ca/notes LOAD \
""
simundump kpool /kinetics/GABA/Ca_ext 0 0.0 0.080003 0.080003 3.1649 3.1649 0 \
  0 39.56 4 /kinetics/geometry 52 yellow 11 3 0
simundump text /kinetics/GABA/Ca_ext/notes 0 ""
call /kinetics/GABA/Ca_ext/notes LOAD \
""
simundump kpool /kinetics/GABA/RR_pool 0 0.0 3 3 118.68 118.68 0 0 39.56 0 \
  /kinetics/geometry 1 yellow 5 1 0
simundump text /kinetics/GABA/RR_pool/notes 0 ""
call /kinetics/GABA/RR_pool/notes LOAD \
""
simundump kpool /kinetics/GABA/vesicle_pool 0 0.0 3.0485 3.0485 120.6 120.6 0 \
  0 39.56 4 /kinetics/geometry 27 yellow 11 -1 0
simundump text /kinetics/GABA/vesicle_pool/notes 0 ""
call /kinetics/GABA/vesicle_pool/notes LOAD \
""
simundump kpool /kinetics/GABA/Docked 0 0.0 0 0 0 0 0 0 39.56 0 \
  /kinetics/geometry 54 blue 9 1 0
simundump text /kinetics/GABA/Docked/notes 0 ""
call /kinetics/GABA/Docked/notes LOAD \
""
simundump kpool /kinetics/GABA/Ca_RR 0 0.0 0 0 0 0 0 0 39.56 0 \
  /kinetics/geometry 51 blue 7 1 0
simundump text /kinetics/GABA/Ca_RR/notes 0 ""
call /kinetics/GABA/Ca_RR/notes LOAD \
""
simundump kpool /kinetics/GABA/GABA 0 0.0 0 0 0 0 0 0 39.56 0 \
  /kinetics/geometry 7 yellow 11 1 0
simundump text /kinetics/GABA/GABA/notes 0 ""
call /kinetics/GABA/GABA/notes LOAD \
""
simundump kreac /kinetics/GABA/remove_Ca 0 0.37917 0.12639 "" white yellow 9 \
  4 0
simundump text /kinetics/GABA/remove_Ca/notes 0 ""
call /kinetics/GABA/remove_Ca/notes LOAD \
""
simundump kreac /kinetics/GABA/remove 0 50000 0 "" white yellow 9 0 0
simundump text /kinetics/GABA/remove/notes 0 ""
call /kinetics/GABA/remove/notes LOAD \
""
simundump kreac /kinetics/GABA/replenish_vesicle 0 1 1 "" white yellow 9 -2 0
simundump text /kinetics/GABA/replenish_vesicle/notes 0 ""
call /kinetics/GABA/replenish_vesicle/notes LOAD \
""
simundump kreac /kinetics/GABA/vesicle_release 0 0.00063898 0 "" white yellow \
  10 2 0
simundump text /kinetics/GABA/vesicle_release/notes 0 \
  "High cooperativity, 4 or higher. Several refs: McDargh and O-Shaughnessy, BioRxiv 2021 Voleti, Jaczynska, Rizo, eLife 2020 Chen.... Scheller, Cell 1999"
call /kinetics/GABA/vesicle_release/notes LOAD \
"High cooperativity, 4 or higher. Several refs: McDargh and O-Shaughnessy, BioRxiv 2021 Voleti, Jaczynska, Rizo, eLife 2020 Chen.... Scheller, Cell 1999"
simundump kreac /kinetics/GABA/Ca_bind_RR 0 0.0063898 2 "" white blue 6 2 0
simundump text /kinetics/GABA/Ca_bind_RR/notes 0 ""
call /kinetics/GABA/Ca_bind_RR/notes LOAD \
""
simundump kreac /kinetics/GABA/docking 0 500 0 "" white blue 8 2 0
simundump text /kinetics/GABA/docking/notes 0 ""
call /kinetics/GABA/docking/notes LOAD \
""
simundump kreac /kinetics/GABA/undocking 0 500 0 "" white blue 7 0 0
simundump text /kinetics/GABA/undocking/notes 0 ""
call /kinetics/GABA/undocking/notes LOAD \
""
simundump kpool /kinetics/GABA/Buffer 0 0 5 5 197.07 197.07 0 0 39.415 0 \
  /kinetics/geometry[1] 47 blue 3 3 0
simundump text /kinetics/GABA/Buffer/notes 0 ""
call /kinetics/GABA/Buffer/notes LOAD \
""
simundump kpool /kinetics/GABA/Ca4.Buffer 0 0 0 0 0 0 0 0 39.415 0 \
  /kinetics/geometry[1] 50 blue 3 1 0
simundump text /kinetics/GABA/Ca4.Buffer/notes 0 ""
call /kinetics/GABA/Ca4.Buffer/notes LOAD \
""
simundump kreac /kinetics/GABA/Ca_bind_buffer 0 1.6573e-05 1000 "" white blue \
  4 2 0
simundump text /kinetics/GABA/Ca_bind_buffer/notes 0 ""
call /kinetics/GABA/Ca_bind_buffer/notes LOAD \
""
simundump text /kinetics/geometry/notes 0 ""
call /kinetics/geometry/notes LOAD \
""
simundump text /kinetics/geometry[1]/notes 0 ""
call /kinetics/geometry[1]/notes LOAD \
""
simundump xgraph /graphs/conc1 0 0 99 0.001 0.999 0
simundump xgraph /graphs/conc2 0 0 100 0 1 0
simundump xplot /graphs/conc1/GABA.Co 3 524288 \
  "delete_plot.w <s> <d>; edit_plot.D <w>" 7 0 0 1
simundump xplot /graphs/conc2/glu.Co 3 524288 \
  "delete_plot.w <s> <d>; edit_plot.D <w>" 7 0 0 1
simundump xgraph /moregraphs/conc3 0 0 100 0 1 0
simundump xgraph /moregraphs/conc4 0 0 100 0 1 0
simundump xcoredraw /edit/draw 0 1 26 -4 7
simundump xtree /edit/draw/tree 0 \
  /kinetics/#[],/kinetics/#[]/#[],/kinetics/#[]/#[]/#[][TYPE!=proto],/kinetics/#[]/#[]/#[][TYPE!=linkinfo]/##[] \
  "edit_elm.D <v>; drag_from_edit.w <d> <S> <x> <y> <z>" auto 0.6
simundump xtext /file/notes 0 1
xtextload /file/notes \
"Based on BothPresyn77.g" \
"27 Apr 20204"
addmsg /kinetics/glu/remove_Ca /kinetics/glu/Ca REAC A B 
addmsg /kinetics/glu/remove_Ca /kinetics/glu/Ca REAC A B 
addmsg /kinetics/glu/remove /kinetics/glu/Ca REAC B A 
addmsg /kinetics/glu/remove /kinetics/glu/Ca REAC B A 
addmsg /kinetics/glu/vesicle_release /kinetics/glu/Ca REAC A B 
addmsg /kinetics/glu/vesicle_release /kinetics/glu/Ca REAC A B 
addmsg /kinetics/glu/Ca_bind_RR /kinetics/glu/Ca REAC A B 
addmsg /kinetics/glu/Ca_bind_RR /kinetics/glu/Ca REAC A B 
addmsg /kinetics/glu/docking /kinetics/glu/Ca REAC B A 
addmsg /kinetics/glu/docking /kinetics/glu/Ca REAC B A 
addmsg /kinetics/glu/Ca_bind_buffer /kinetics/glu/Ca REAC A B 
addmsg /kinetics/glu/Ca_bind_buffer /kinetics/glu/Ca REAC A B 
addmsg /kinetics/glu/Ca_bind_buffer /kinetics/glu/Ca REAC A B 
addmsg /kinetics/glu/Ca_bind_buffer /kinetics/glu/Ca REAC A B 
addmsg /kinetics/glu/remove_Ca /kinetics/glu/Ca_ext REAC B A 
addmsg /kinetics/glu/remove_Ca /kinetics/glu/Ca_ext REAC B A 
addmsg /kinetics/glu/replenish_vesicle /kinetics/glu/RR_pool REAC B A 
addmsg /kinetics/glu/Ca_bind_RR /kinetics/glu/RR_pool REAC A B 
addmsg /kinetics/glu/undocking /kinetics/glu/RR_pool REAC B A 
addmsg /kinetics/glu/remove /kinetics/glu/vesicle_pool REAC B A 
addmsg /kinetics/glu/replenish_vesicle /kinetics/glu/vesicle_pool REAC A B 
addmsg /kinetics/glu/remove /kinetics/glu/glu REAC A B 
addmsg /kinetics/glu/vesicle_release /kinetics/glu/glu REAC B A 
addmsg /kinetics/glu/vesicle_release /kinetics/glu/Docked REAC A B 
addmsg /kinetics/glu/docking /kinetics/glu/Docked REAC B A 
addmsg /kinetics/glu/undocking /kinetics/glu/Docked REAC A B 
addmsg /kinetics/glu/Ca_bind_RR /kinetics/glu/Ca_RR REAC B A 
addmsg /kinetics/glu/docking /kinetics/glu/Ca_RR REAC A B 
addmsg /kinetics/glu/Ca /kinetics/glu/remove_Ca SUBSTRATE n 
addmsg /kinetics/glu/Ca /kinetics/glu/remove_Ca SUBSTRATE n 
addmsg /kinetics/glu/Ca_ext /kinetics/glu/remove_Ca PRODUCT n 
addmsg /kinetics/glu/Ca_ext /kinetics/glu/remove_Ca PRODUCT n 
addmsg /kinetics/glu/glu /kinetics/glu/remove SUBSTRATE n 
addmsg /kinetics/glu/Ca /kinetics/glu/remove PRODUCT n 
addmsg /kinetics/glu/Ca /kinetics/glu/remove PRODUCT n 
addmsg /kinetics/glu/vesicle_pool /kinetics/glu/remove PRODUCT n 
addmsg /kinetics/glu/vesicle_pool /kinetics/glu/replenish_vesicle SUBSTRATE n 
addmsg /kinetics/glu/RR_pool /kinetics/glu/replenish_vesicle PRODUCT n 
addmsg /kinetics/glu/Ca /kinetics/glu/vesicle_release SUBSTRATE n 
addmsg /kinetics/glu/Ca /kinetics/glu/vesicle_release SUBSTRATE n 
addmsg /kinetics/glu/Docked /kinetics/glu/vesicle_release SUBSTRATE n 
addmsg /kinetics/glu/glu /kinetics/glu/vesicle_release PRODUCT n 
addmsg /kinetics/glu/Ca /kinetics/glu/Ca_bind_RR SUBSTRATE n 
addmsg /kinetics/glu/Ca /kinetics/glu/Ca_bind_RR SUBSTRATE n 
addmsg /kinetics/glu/RR_pool /kinetics/glu/Ca_bind_RR SUBSTRATE n 
addmsg /kinetics/glu/Ca_RR /kinetics/glu/Ca_bind_RR PRODUCT n 
addmsg /kinetics/glu/Ca_RR /kinetics/glu/docking SUBSTRATE n 
addmsg /kinetics/glu/Ca /kinetics/glu/docking PRODUCT n 
addmsg /kinetics/glu/Ca /kinetics/glu/docking PRODUCT n 
addmsg /kinetics/glu/Docked /kinetics/glu/docking PRODUCT n 
addmsg /kinetics/glu/Docked /kinetics/glu/undocking SUBSTRATE n 
addmsg /kinetics/glu/RR_pool /kinetics/glu/undocking PRODUCT n 
addmsg /kinetics/glu/Ca_bind_buffer /kinetics/glu/Buffer REAC A B 
addmsg /kinetics/glu/Ca_bind_buffer /kinetics/glu/Ca4.Buffer REAC B A 
addmsg /kinetics/glu/Buffer /kinetics/glu/Ca_bind_buffer SUBSTRATE n 
addmsg /kinetics/glu/Ca /kinetics/glu/Ca_bind_buffer SUBSTRATE n 
addmsg /kinetics/glu/Ca /kinetics/glu/Ca_bind_buffer SUBSTRATE n 
addmsg /kinetics/glu/Ca /kinetics/glu/Ca_bind_buffer SUBSTRATE n 
addmsg /kinetics/glu/Ca /kinetics/glu/Ca_bind_buffer SUBSTRATE n 
addmsg /kinetics/glu/Ca4.Buffer /kinetics/glu/Ca_bind_buffer PRODUCT n 
addmsg /kinetics/GABA/remove_Ca /kinetics/GABA/Ca REAC A B 
addmsg /kinetics/GABA/remove_Ca /kinetics/GABA/Ca REAC A B 
addmsg /kinetics/GABA/remove /kinetics/GABA/Ca REAC B A 
addmsg /kinetics/GABA/remove /kinetics/GABA/Ca REAC B A 
addmsg /kinetics/GABA/vesicle_release /kinetics/GABA/Ca REAC A B 
addmsg /kinetics/GABA/vesicle_release /kinetics/GABA/Ca REAC A B 
addmsg /kinetics/GABA/Ca_bind_RR /kinetics/GABA/Ca REAC A B 
addmsg /kinetics/GABA/Ca_bind_RR /kinetics/GABA/Ca REAC A B 
addmsg /kinetics/GABA/docking /kinetics/GABA/Ca REAC B A 
addmsg /kinetics/GABA/docking /kinetics/GABA/Ca REAC B A 
addmsg /kinetics/GABA/Ca_bind_buffer /kinetics/GABA/Ca REAC A B 
addmsg /kinetics/GABA/Ca_bind_buffer /kinetics/GABA/Ca REAC A B 
addmsg /kinetics/GABA/Ca_bind_buffer /kinetics/GABA/Ca REAC A B 
addmsg /kinetics/GABA/Ca_bind_buffer /kinetics/GABA/Ca REAC A B 
addmsg /kinetics/GABA/remove_Ca /kinetics/GABA/Ca_ext REAC B A 
addmsg /kinetics/GABA/remove_Ca /kinetics/GABA/Ca_ext REAC B A 
addmsg /kinetics/GABA/replenish_vesicle /kinetics/GABA/RR_pool REAC B A 
addmsg /kinetics/GABA/Ca_bind_RR /kinetics/GABA/RR_pool REAC A B 
addmsg /kinetics/GABA/undocking /kinetics/GABA/RR_pool REAC B A 
addmsg /kinetics/GABA/remove /kinetics/GABA/vesicle_pool REAC B A 
addmsg /kinetics/GABA/replenish_vesicle /kinetics/GABA/vesicle_pool REAC A B 
addmsg /kinetics/GABA/vesicle_release /kinetics/GABA/Docked REAC A B 
addmsg /kinetics/GABA/docking /kinetics/GABA/Docked REAC B A 
addmsg /kinetics/GABA/undocking /kinetics/GABA/Docked REAC A B 
addmsg /kinetics/GABA/Ca_bind_RR /kinetics/GABA/Ca_RR REAC B A 
addmsg /kinetics/GABA/docking /kinetics/GABA/Ca_RR REAC A B 
addmsg /kinetics/GABA/remove /kinetics/GABA/GABA REAC A B 
addmsg /kinetics/GABA/vesicle_release /kinetics/GABA/GABA REAC B A 
addmsg /kinetics/GABA/Ca /kinetics/GABA/remove_Ca SUBSTRATE n 
addmsg /kinetics/GABA/Ca /kinetics/GABA/remove_Ca SUBSTRATE n 
addmsg /kinetics/GABA/Ca_ext /kinetics/GABA/remove_Ca PRODUCT n 
addmsg /kinetics/GABA/Ca_ext /kinetics/GABA/remove_Ca PRODUCT n 
addmsg /kinetics/GABA/GABA /kinetics/GABA/remove SUBSTRATE n 
addmsg /kinetics/GABA/Ca /kinetics/GABA/remove PRODUCT n 
addmsg /kinetics/GABA/Ca /kinetics/GABA/remove PRODUCT n 
addmsg /kinetics/GABA/vesicle_pool /kinetics/GABA/remove PRODUCT n 
addmsg /kinetics/GABA/vesicle_pool /kinetics/GABA/replenish_vesicle SUBSTRATE n 
addmsg /kinetics/GABA/RR_pool /kinetics/GABA/replenish_vesicle PRODUCT n 
addmsg /kinetics/GABA/Ca /kinetics/GABA/vesicle_release SUBSTRATE n 
addmsg /kinetics/GABA/Ca /kinetics/GABA/vesicle_release SUBSTRATE n 
addmsg /kinetics/GABA/Docked /kinetics/GABA/vesicle_release SUBSTRATE n 
addmsg /kinetics/GABA/GABA /kinetics/GABA/vesicle_release PRODUCT n 
addmsg /kinetics/GABA/Ca /kinetics/GABA/Ca_bind_RR SUBSTRATE n 
addmsg /kinetics/GABA/Ca /kinetics/GABA/Ca_bind_RR SUBSTRATE n 
addmsg /kinetics/GABA/RR_pool /kinetics/GABA/Ca_bind_RR SUBSTRATE n 
addmsg /kinetics/GABA/Ca_RR /kinetics/GABA/Ca_bind_RR PRODUCT n 
addmsg /kinetics/GABA/Ca_RR /kinetics/GABA/docking SUBSTRATE n 
addmsg /kinetics/GABA/Ca /kinetics/GABA/docking PRODUCT n 
addmsg /kinetics/GABA/Ca /kinetics/GABA/docking PRODUCT n 
addmsg /kinetics/GABA/Docked /kinetics/GABA/docking PRODUCT n 
addmsg /kinetics/GABA/Docked /kinetics/GABA/undocking SUBSTRATE n 
addmsg /kinetics/GABA/RR_pool /kinetics/GABA/undocking PRODUCT n 
addmsg /kinetics/GABA/Ca_bind_buffer /kinetics/GABA/Buffer REAC A B 
addmsg /kinetics/GABA/Ca_bind_buffer /kinetics/GABA/Ca4.Buffer REAC B A 
addmsg /kinetics/GABA/Buffer /kinetics/GABA/Ca_bind_buffer SUBSTRATE n 
addmsg /kinetics/GABA/Ca /kinetics/GABA/Ca_bind_buffer SUBSTRATE n 
addmsg /kinetics/GABA/Ca /kinetics/GABA/Ca_bind_buffer SUBSTRATE n 
addmsg /kinetics/GABA/Ca /kinetics/GABA/Ca_bind_buffer SUBSTRATE n 
addmsg /kinetics/GABA/Ca /kinetics/GABA/Ca_bind_buffer SUBSTRATE n 
addmsg /kinetics/GABA/Ca4.Buffer /kinetics/GABA/Ca_bind_buffer PRODUCT n 
addmsg /kinetics/GABA/GABA /graphs/conc1/GABA.Co PLOT Co *GABA.Co *7 
addmsg /kinetics/glu/glu /graphs/conc2/glu.Co PLOT Co *glu.Co *7 
enddump
// End of dump

call /kinetics/glu/notes LOAD \
"Compartment"
call /kinetics/glu/Ca/notes LOAD \
"Compartment"
call /kinetics/glu/vesicle_release/notes LOAD \
"High cooperativity, 4 or higher. Several refs: McDargh and O-Shaughnessy, BioRxiv 2021 Voleti, Jaczynska, Rizo, eLife 2020 Chen.... Scheller, Cell 1999"
call /kinetics/GABA/notes LOAD \
"Compartment"
call /kinetics/GABA/vesicle_release/notes LOAD \
"High cooperativity, 4 or higher. Several refs: McDargh and O-Shaughnessy, BioRxiv 2021 Voleti, Jaczynska, Rizo, eLife 2020 Chen.... Scheller, Cell 1999"
complete_loading
