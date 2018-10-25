# Compile postAnalyzers

./rootcom postAnalyzer_ZnnG_qcd_above0p5 ZnnG_qcd_above0p5
./rootcom postAnalyzer_ZnnG_qcd_below0p5 ZnnG_qcd_below0p5
./rootcom postAnalyzer_ZeeG_qcd ZeeG_qcd
./rootcom postAnalyzer_ZmmG_qcd ZmmG_qcd
./rootcom postAnalyzer_WenG_qcd WenG_qcd
./rootcom postAnalyzer_WmnG_qcd WmnG_qcd


# Run postAnalyzers to get expected yield distributions and their systematic shifts

./MakeCondorFiles.csh ZnnG_qcd_above0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016B_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016B_v1/170714_132232/0000/  ZnnG_qcd_b0000_above0p5.root -1 1000 ZnnG_qcd_b0000_above0p5
./MakeCondorFiles.csh ZnnG_qcd_above0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016B_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016B_v1/170714_132232/0001/  ZnnG_qcd_b0001_above0p5.root -1 1000 ZnnG_qcd_b0001_above0p5
./MakeCondorFiles.csh ZnnG_qcd_above0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016B_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016B_v1/170714_132232/0002/  ZnnG_qcd_b0002_above0p5.root -1 1000 ZnnG_qcd_b0002_above0p5
./MakeCondorFiles.csh ZnnG_qcd_above0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016B_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016B_v1/170714_132232/0003/  ZnnG_qcd_b0003_above0p5.root -1 1000 ZnnG_qcd_b0003_above0p5
./MakeCondorFiles.csh ZnnG_qcd_above0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016C_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016C_v1/170714_132346/0000/  ZnnG_qcd_c0000_above0p5.root -1 1000 ZnnG_qcd_c0000_above0p5
./MakeCondorFiles.csh ZnnG_qcd_above0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016C_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016C_v1/170714_132346/0001/  ZnnG_qcd_c0001_above0p5.root -1 1000 ZnnG_qcd_c0001_above0p5
./MakeCondorFiles.csh ZnnG_qcd_above0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016D_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016D_v1/170714_132504/0000/  ZnnG_qcd_d0000_above0p5.root -1 1000 ZnnG_qcd_d0000_above0p5
./MakeCondorFiles.csh ZnnG_qcd_above0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016D_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016D_v1/170714_132504/0001/  ZnnG_qcd_d0001_above0p5.root -1 1000 ZnnG_qcd_d0001_above0p5
./MakeCondorFiles.csh ZnnG_qcd_above0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016E_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016E_v1/170714_160205/0000/  ZnnG_qcd_e0000_above0p5.root -1 1000 ZnnG_qcd_e0000_above0p5
./MakeCondorFiles.csh ZnnG_qcd_above0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016E_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016E_v1/170714_160205/0001/  ZnnG_qcd_e0001_above0p5.root -1 1000 ZnnG_qcd_e0001_above0p5
./MakeCondorFiles.csh ZnnG_qcd_above0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016F_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016F_v1/170714_160448/0000/  ZnnG_qcd_f0000_above0p5.root -1 1000 ZnnG_qcd_f0000_above0p5
./MakeCondorFiles.csh ZnnG_qcd_above0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016F_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016F_v1/170714_160448/0001/  ZnnG_qcd_f0001_above0p5.root -1 1000 ZnnG_qcd_f0001_above0p5
./MakeCondorFiles.csh ZnnG_qcd_above0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016G_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016G_v1/170714_160541/0000/  ZnnG_qcd_g0000_above0p5.root -1 1000 ZnnG_qcd_g0000_above0p5
./MakeCondorFiles.csh ZnnG_qcd_above0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016G_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016G_v1/170714_160541/0001/  ZnnG_qcd_g0001_above0p5.root -1 1000 ZnnG_qcd_g0001_above0p5
./MakeCondorFiles.csh ZnnG_qcd_above0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016G_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016G_v1/170714_160541/0002/  ZnnG_qcd_g0002_above0p5.root -1 1000 ZnnG_qcd_g0002_above0p5
./MakeCondorFiles.csh ZnnG_qcd_above0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver2_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver2_v1/170714_160802/0000/  ZnnG_qcd_h0000_above0p5.root -1 1000 ZnnG_qcd_h0000_above0p5
./MakeCondorFiles.csh ZnnG_qcd_above0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver2_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver2_v1/170714_160802/0001/  ZnnG_qcd_h0001_above0p5.root -1 1000 ZnnG_qcd_h0001_above0p5
./MakeCondorFiles.csh ZnnG_qcd_above0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver2_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver2_v1/170714_160802/0002/  ZnnG_qcd_h0002_above0p5.root -1 1000 ZnnG_qcd_h0002_above0p5
./MakeCondorFiles.csh ZnnG_qcd_above0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver2_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver2_v1/170714_160802/0003/  ZnnG_qcd_h0003_above0p5.root -1 1000 ZnnG_qcd_h0003_above0p5
./MakeCondorFiles.csh ZnnG_qcd_above0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver3_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver3_v1/170714_160847/0000/  ZnnG_qcd_h0004_above0p5.root -1 1000 ZnnG_qcd_h0004_above0p5

./MakeCondorFiles.csh ZnnG_qcd_below0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016B_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016B_v1/170714_132232/0000/  ZnnG_qcd_b0000_below0p5.root -1 1000 ZnnG_qcd_b0000_below0p5
./MakeCondorFiles.csh ZnnG_qcd_below0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016B_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016B_v1/170714_132232/0001/  ZnnG_qcd_b0001_below0p5.root -1 1000 ZnnG_qcd_b0001_below0p5
./MakeCondorFiles.csh ZnnG_qcd_below0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016B_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016B_v1/170714_132232/0002/  ZnnG_qcd_b0002_below0p5.root -1 1000 ZnnG_qcd_b0002_below0p5
./MakeCondorFiles.csh ZnnG_qcd_below0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016B_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016B_v1/170714_132232/0003/  ZnnG_qcd_b0003_below0p5.root -1 1000 ZnnG_qcd_b0003_below0p5
./MakeCondorFiles.csh ZnnG_qcd_below0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016C_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016C_v1/170714_132346/0000/  ZnnG_qcd_c0000_below0p5.root -1 1000 ZnnG_qcd_c0000_below0p5
./MakeCondorFiles.csh ZnnG_qcd_below0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016C_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016C_v1/170714_132346/0001/  ZnnG_qcd_c0001_below0p5.root -1 1000 ZnnG_qcd_c0001_below0p5
./MakeCondorFiles.csh ZnnG_qcd_below0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016D_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016D_v1/170714_132504/0000/  ZnnG_qcd_d0000_below0p5.root -1 1000 ZnnG_qcd_d0000_below0p5
./MakeCondorFiles.csh ZnnG_qcd_below0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016D_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016D_v1/170714_132504/0001/  ZnnG_qcd_d0001_below0p5.root -1 1000 ZnnG_qcd_d0001_below0p5
./MakeCondorFiles.csh ZnnG_qcd_below0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016E_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016E_v1/170714_160205/0000/  ZnnG_qcd_e0000_below0p5.root -1 1000 ZnnG_qcd_e0000_below0p5
./MakeCondorFiles.csh ZnnG_qcd_below0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016E_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016E_v1/170714_160205/0001/  ZnnG_qcd_e0001_below0p5.root -1 1000 ZnnG_qcd_e0001_below0p5
./MakeCondorFiles.csh ZnnG_qcd_below0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016F_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016F_v1/170714_160448/0000/  ZnnG_qcd_f0000_below0p5.root -1 1000 ZnnG_qcd_f0000_below0p5
./MakeCondorFiles.csh ZnnG_qcd_below0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016F_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016F_v1/170714_160448/0001/  ZnnG_qcd_f0001_below0p5.root -1 1000 ZnnG_qcd_f0001_below0p5
./MakeCondorFiles.csh ZnnG_qcd_below0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016G_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016G_v1/170714_160541/0000/  ZnnG_qcd_g0000_below0p5.root -1 1000 ZnnG_qcd_g0000_below0p5
./MakeCondorFiles.csh ZnnG_qcd_below0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016G_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016G_v1/170714_160541/0001/  ZnnG_qcd_g0001_below0p5.root -1 1000 ZnnG_qcd_g0001_below0p5
./MakeCondorFiles.csh ZnnG_qcd_below0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016G_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016G_v1/170714_160541/0002/  ZnnG_qcd_g0002_below0p5.root -1 1000 ZnnG_qcd_g0002_below0p5
./MakeCondorFiles.csh ZnnG_qcd_below0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver2_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver2_v1/170714_160802/0000/  ZnnG_qcd_h0000_below0p5.root -1 1000 ZnnG_qcd_h0000_below0p5
./MakeCondorFiles.csh ZnnG_qcd_below0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver2_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver2_v1/170714_160802/0001/  ZnnG_qcd_h0001_below0p5.root -1 1000 ZnnG_qcd_h0001_below0p5
./MakeCondorFiles.csh ZnnG_qcd_below0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver2_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver2_v1/170714_160802/0002/  ZnnG_qcd_h0002_below0p5.root -1 1000 ZnnG_qcd_h0002_below0p5
./MakeCondorFiles.csh ZnnG_qcd_below0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver2_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver2_v1/170714_160802/0003/  ZnnG_qcd_h0003_below0p5.root -1 1000 ZnnG_qcd_h0003_below0p5
./MakeCondorFiles.csh ZnnG_qcd_below0p5 /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver3_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver3_v1/170714_160847/0000/  ZnnG_qcd_h0004_below0p5.root -1 1000 ZnnG_qcd_h0004_below0p5

./MakeCondorFiles.csh ZeeG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016B_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016B_v1/170714_132232/0000/  ZeeG_qcd_b0000.root -1 1000 ZeeG_qcd_b0000
./MakeCondorFiles.csh ZeeG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016B_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016B_v1/170714_132232/0001/  ZeeG_qcd_b0001.root -1 1000 ZeeG_qcd_b0001
./MakeCondorFiles.csh ZeeG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016B_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016B_v1/170714_132232/0002/  ZeeG_qcd_b0002.root -1 1000 ZeeG_qcd_b0002
./MakeCondorFiles.csh ZeeG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016B_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016B_v1/170714_132232/0003/  ZeeG_qcd_b0003.root -1 1000 ZeeG_qcd_b0003
./MakeCondorFiles.csh ZeeG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016C_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016C_v1/170714_132346/0000/  ZeeG_qcd_c0000.root -1 1000 ZeeG_qcd_c0000
./MakeCondorFiles.csh ZeeG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016C_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016C_v1/170714_132346/0001/  ZeeG_qcd_c0001.root -1 1000 ZeeG_qcd_c0001
./MakeCondorFiles.csh ZeeG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016D_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016D_v1/170714_132504/0000/  ZeeG_qcd_d0000.root -1 1000 ZeeG_qcd_d0000
./MakeCondorFiles.csh ZeeG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016D_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016D_v1/170714_132504/0001/  ZeeG_qcd_d0001.root -1 1000 ZeeG_qcd_d0001
./MakeCondorFiles.csh ZeeG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016E_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016E_v1/170714_160205/0000/  ZeeG_qcd_e0000.root -1 1000 ZeeG_qcd_e0000
./MakeCondorFiles.csh ZeeG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016E_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016E_v1/170714_160205/0001/  ZeeG_qcd_e0001.root -1 1000 ZeeG_qcd_e0001
./MakeCondorFiles.csh ZeeG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016F_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016F_v1/170714_160448/0000/  ZeeG_qcd_f0000.root -1 1000 ZeeG_qcd_f0000
./MakeCondorFiles.csh ZeeG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016F_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016F_v1/170714_160448/0001/  ZeeG_qcd_f0001.root -1 1000 ZeeG_qcd_f0001
./MakeCondorFiles.csh ZeeG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016G_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016G_v1/170714_160541/0000/  ZeeG_qcd_g0000.root -1 1000 ZeeG_qcd_g0000
./MakeCondorFiles.csh ZeeG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016G_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016G_v1/170714_160541/0001/  ZeeG_qcd_g0001.root -1 1000 ZeeG_qcd_g0001
./MakeCondorFiles.csh ZeeG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016G_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016G_v1/170714_160541/0002/  ZeeG_qcd_g0002.root -1 1000 ZeeG_qcd_g0002
./MakeCondorFiles.csh ZeeG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver2_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver2_v1/170714_160802/0000/  ZeeG_qcd_h0000.root -1 1000 ZeeG_qcd_h0000
./MakeCondorFiles.csh ZeeG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver2_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver2_v1/170714_160802/0001/  ZeeG_qcd_h0001.root -1 1000 ZeeG_qcd_h0001
./MakeCondorFiles.csh ZeeG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver2_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver2_v1/170714_160802/0002/  ZeeG_qcd_h0002.root -1 1000 ZeeG_qcd_h0002
./MakeCondorFiles.csh ZeeG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver2_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver2_v1/170714_160802/0003/  ZeeG_qcd_h0003.root -1 1000 ZeeG_qcd_h0003
./MakeCondorFiles.csh ZeeG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver3_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver3_v1/170714_160847/0000/  ZeeG_qcd_h0004.root -1 1000 ZeeG_qcd_h0004

./MakeCondorFiles.csh ZmmG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016B_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016B_v1/170714_132232/0000/  ZmmG_qcd_b0000.root -1 1000 ZmmG_qcd_b0000
./MakeCondorFiles.csh ZmmG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016B_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016B_v1/170714_132232/0001/  ZmmG_qcd_b0001.root -1 1000 ZmmG_qcd_b0001
./MakeCondorFiles.csh ZmmG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016B_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016B_v1/170714_132232/0002/  ZmmG_qcd_b0002.root -1 1000 ZmmG_qcd_b0002
./MakeCondorFiles.csh ZmmG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016B_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016B_v1/170714_132232/0003/  ZmmG_qcd_b0003.root -1 1000 ZmmG_qcd_b0003
./MakeCondorFiles.csh ZmmG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016C_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016C_v1/170714_132346/0000/  ZmmG_qcd_c0000.root -1 1000 ZmmG_qcd_c0000
./MakeCondorFiles.csh ZmmG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016C_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016C_v1/170714_132346/0001/  ZmmG_qcd_c0001.root -1 1000 ZmmG_qcd_c0001
./MakeCondorFiles.csh ZmmG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016D_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016D_v1/170714_132504/0000/  ZmmG_qcd_d0000.root -1 1000 ZmmG_qcd_d0000
./MakeCondorFiles.csh ZmmG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016D_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016D_v1/170714_132504/0001/  ZmmG_qcd_d0001.root -1 1000 ZmmG_qcd_d0001
./MakeCondorFiles.csh ZmmG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016E_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016E_v1/170714_160205/0000/  ZmmG_qcd_e0000.root -1 1000 ZmmG_qcd_e0000
./MakeCondorFiles.csh ZmmG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016E_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016E_v1/170714_160205/0001/  ZmmG_qcd_e0001.root -1 1000 ZmmG_qcd_e0001
./MakeCondorFiles.csh ZmmG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016F_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016F_v1/170714_160448/0000/  ZmmG_qcd_f0000.root -1 1000 ZmmG_qcd_f0000
./MakeCondorFiles.csh ZmmG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016F_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016F_v1/170714_160448/0001/  ZmmG_qcd_f0001.root -1 1000 ZmmG_qcd_f0001
./MakeCondorFiles.csh ZmmG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016G_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016G_v1/170714_160541/0000/  ZmmG_qcd_g0000.root -1 1000 ZmmG_qcd_g0000
./MakeCondorFiles.csh ZmmG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016G_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016G_v1/170714_160541/0001/  ZmmG_qcd_g0001.root -1 1000 ZmmG_qcd_g0001
./MakeCondorFiles.csh ZmmG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016G_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016G_v1/170714_160541/0002/  ZmmG_qcd_g0002.root -1 1000 ZmmG_qcd_g0002
./MakeCondorFiles.csh ZmmG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver2_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver2_v1/170714_160802/0000/  ZmmG_qcd_h0000.root -1 1000 ZmmG_qcd_h0000
./MakeCondorFiles.csh ZmmG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver2_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver2_v1/170714_160802/0001/  ZmmG_qcd_h0001.root -1 1000 ZmmG_qcd_h0001
./MakeCondorFiles.csh ZmmG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver2_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver2_v1/170714_160802/0002/  ZmmG_qcd_h0002.root -1 1000 ZmmG_qcd_h0002
./MakeCondorFiles.csh ZmmG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver2_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver2_v1/170714_160802/0003/  ZmmG_qcd_h0003.root -1 1000 ZmmG_qcd_h0003
./MakeCondorFiles.csh ZmmG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver3_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver3_v1/170714_160847/0000/  ZmmG_qcd_h0004.root -1 1000 ZmmG_qcd_h0004

./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016B_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016B_v1/170714_132232/0000/  WenG_qcd_b0000.root -1 1000 WenG_qcd_b0000
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016B_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016B_v1/170714_132232/0001/  WenG_qcd_b0001.root -1 1000 WenG_qcd_b0001
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016B_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016B_v1/170714_132232/0002/  WenG_qcd_b0002.root -1 1000 WenG_qcd_b0002
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016B_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016B_v1/170714_132232/0003/  WenG_qcd_b0003.root -1 1000 WenG_qcd_b0003
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016C_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016C_v1/170714_132346/0000/  WenG_qcd_c0000.root -1 1000 WenG_qcd_c0000
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016C_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016C_v1/170714_132346/0001/  WenG_qcd_c0001.root -1 1000 WenG_qcd_c0001
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016D_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016D_v1/170714_132504/0000/  WenG_qcd_d0000.root -1 1000 WenG_qcd_d0000
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016D_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016D_v1/170714_132504/0001/  WenG_qcd_d0001.root -1 1000 WenG_qcd_d0001
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016E_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016E_v1/170714_160205/0000/  WenG_qcd_e0000.root -1 1000 WenG_qcd_e0000
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016E_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016E_v1/170714_160205/0001/  WenG_qcd_e0001.root -1 1000 WenG_qcd_e0001
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016F_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016F_v1/170714_160448/0000/  WenG_qcd_f0000.root -1 1000 WenG_qcd_f0000
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016F_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016F_v1/170714_160448/0001/  WenG_qcd_f0001.root -1 1000 WenG_qcd_f0001
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016G_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016G_v1/170714_160541/0000/  WenG_qcd_g0000.root -1 1000 WenG_qcd_g0000
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016G_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016G_v1/170714_160541/0001/  WenG_qcd_g0001.root -1 1000 WenG_qcd_g0001
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016G_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016G_v1/170714_160541/0002/  WenG_qcd_g0002.root -1 1000 WenG_qcd_g0002
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver2_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver2_v1/170714_160802/0000/  WenG_qcd_h0000.root -1 1000 WenG_qcd_h0000
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver2_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver2_v1/170714_160802/0001/  WenG_qcd_h0001.root -1 1000 WenG_qcd_h0001
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver2_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver2_v1/170714_160802/0002/  WenG_qcd_h0002.root -1 1000 WenG_qcd_h0002
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver2_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver2_v1/170714_160802/0003/  WenG_qcd_h0003.root -1 1000 WenG_qcd_h0003
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver3_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver3_v1/170714_160847/0000/  WenG_qcd_h0004.root -1 1000 WenG_qcd_h0004

./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016B_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016B_v1/170714_132232/0000/  WmnG_qcd_b0000.root -1 1000 WmnG_qcd_b0000
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016B_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016B_v1/170714_132232/0001/  WmnG_qcd_b0001.root -1 1000 WmnG_qcd_b0001
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016B_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016B_v1/170714_132232/0002/  WmnG_qcd_b0002.root -1 1000 WmnG_qcd_b0002
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016B_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016B_v1/170714_132232/0003/  WmnG_qcd_b0003.root -1 1000 WmnG_qcd_b0003
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016C_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016C_v1/170714_132346/0000/  WmnG_qcd_c0000.root -1 1000 WmnG_qcd_c0000
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016C_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016C_v1/170714_132346/0001/  WmnG_qcd_c0001.root -1 1000 WmnG_qcd_c0001
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016D_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016D_v1/170714_132504/0000/  WmnG_qcd_d0000.root -1 1000 WmnG_qcd_d0000
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016D_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016D_v1/170714_132504/0001/  WmnG_qcd_d0001.root -1 1000 WmnG_qcd_d0001
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016E_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016E_v1/170714_160205/0000/  WmnG_qcd_e0000.root -1 1000 WmnG_qcd_e0000
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016E_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016E_v1/170714_160205/0001/  WmnG_qcd_e0001.root -1 1000 WmnG_qcd_e0001
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016F_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016F_v1/170714_160448/0000/  WmnG_qcd_f0000.root -1 1000 WmnG_qcd_f0000
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016F_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016F_v1/170714_160448/0001/  WmnG_qcd_f0001.root -1 1000 WmnG_qcd_f0001
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016G_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016G_v1/170714_160541/0000/  WmnG_qcd_g0000.root -1 1000 WmnG_qcd_g0000
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016G_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016G_v1/170714_160541/0001/  WmnG_qcd_g0001.root -1 1000 WmnG_qcd_g0001
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016G_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016G_v1/170714_160541/0002/  WmnG_qcd_g0002.root -1 1000 WmnG_qcd_g0002
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver2_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver2_v1/170714_160802/0000/  WmnG_qcd_h0000.root -1 1000 WmnG_qcd_h0000
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver2_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver2_v1/170714_160802/0001/  WmnG_qcd_h0001.root -1 1000 WmnG_qcd_h0001
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver2_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver2_v1/170714_160802/0002/  WmnG_qcd_h0002.root -1 1000 WmnG_qcd_h0002
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver2_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver2_v1/170714_160802/0003/  WmnG_qcd_h0003.root -1 1000 WmnG_qcd_h0003
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/2016ReMiniAODnutples_14July2017/SinglePhoton_2016H_ver3_v1/SinglePhoton/crab_job_Singlephoton_reminiAOD_2016H_ver3_v1/170714_160847/0000/  WmnG_qcd_h0004.root -1 1000 WmnG_qcd_h0004


## Putting it all together ##
hadd -f ZnnG_qcd_above0p5_all.root ZnnG_qcd*000*_above0p5.root
hadd -f ZnnG_qcd_below0p5_all.root ZnnG_qcd*000*_below0p5.root
hadd -f ZeeG_qcd_all.root ZeeG_qcd*000*.root
hadd -f ZmmG_qcd_all.root ZmmG_qcd*000*.root
hadd -f WenG_qcd_all.root WenG_qcd*000*.root
hadd -f WmnG_qcd_all.root WmnG_qcd*000*.root

# Once these are all successfully merged, delete the separate files to save space
rm ZnnG_qcd*000*_above0p5.root
rm ZnnG_qcd*000*_below0p5.root
rm ZeeG_qcd*000*.root
rm ZmmG_qcd*000*.root
rm WenG_qcd*000*.root
rm WmnG_qcd*000*.root


## zeeg_plotter_qcd.C shows how to assemble the final yield histograms and find the systematic uncertainties
## It takes two inputs: ZeeG_qcd_all.root and ZeeG_data_all.root
## If you just want to run a quick test, you could make a temporary ZeeG_data_all.root by copying ZeeG_qcd_all.root:
# cp ZeeG_qcd_all.root ZeeG_data_all.root
root -l -q -b zeeg_plotter_qcd.C

