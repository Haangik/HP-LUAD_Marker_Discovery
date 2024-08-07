### Previous markers were manually curated and organized to character vectors
## CCP markers, Bueno et al. 2015
### KIAA0101 -> PCLAF, ORC6L -> ORC6, CDC2 -> CDK1
CCP.markers<-c("ASPM", "CDCA8", "MCM10", "FOXM1", "CDC20", "CDKN3", "BIRC5", "DLGAP5", "KIF20A", "BUB1B", "PRC1",
               "TK1", "CEP55", "PBK", "RAD54L", "NUSAP1", "RRM2", "PCLAF", "ORC6", "RAD51", "CENPM", 
               "SKA1", "CENPF", "KIF11", "PTTG1", "CDK1", "DTL", "PLK1", "CDCA3", "ASF1B", "TOP2A")
## CIN markers
### CIN-25, Carter et al. 2006
### H2AFZ -> H2AZ1
CIN25<-c('PCNA', 'H2AZ1', 'CCT5', 'MCM7', 'TGIF2', 'NCAPD2', 'FEN1', 'RFC4', 'MCM2', 'UBE2C', 
         'CDK1', 'ESPL1', 'CCNB2', 'TRIP13', 'TPX2', 'TOP2A', 'RAD51AP1', 'MAD2L1', 
         'FOXM1', 'CDC45', 'MELK', 'KIF20A', 'TTK', 'PRC1', 'RNASEH2A')

# ESLA-7, Krzystanek et al. 2016
### FGFR1OP -> CEP43, NCGAP -> NCAPG
ESLA7<-c('ADAM10', 'DLGAP5', 'RAD51AP1', 'CEP43', 'NCAPG', 'KIF15', 'ASPM')

## Zhu et al. 2010
### FAM64A -> PIMREG, IKBKAP -> ELP1
Zhu.15<-c('ATP1B1', 'TRIM14', 'PIMREG', 'FOSL2', 'HEXIM1', 'MB', 'L1CAM', 'UMPS',
          'EDN3', 'STMN2', 'MYT1L', 'ELP1', 'MLANA', 'MDM2', 'ZNF236')

## Park et al. 2012 (193 probes, 145 genes)
Park.2012<-unique(read.table("Park_2012_symbol.txt")$V1)

## Lu et al. 2013
Lu.16<-c('ADD1', 'DNM1L', 'DSG2', 'DSP', 'HMGB2', 'KPNB1', 'LMNB1', 'MAPT', 'OCLN', 'PAK2', 'PKP1', 'PRKCD', 'PRKCQ',
         'SATB1', 'STK24', 'TJP1')

## Okayama et al. 2014
Okayama.4<-c('BRCA1', 'HIF1A', 'DLC1', 'XPO1')

## Gentles et al. 2015
## FAIM3 -> FCMR
Gentles.9<-c('MAD2L1', 'GINS1', 'FCGRT', 'TNIK', 'KDM6A', 'BCAM', 'KRT6A', 'SLC2A1', 'FCMR')

## Li et al. 2018
Li.8<-c('DLGAP5', 'KIF11', 'RAD51AP1', 'CCNB1', 'AURKA', 'CDC6', 'OIP5', 'NCAPG')

## Songyang et al. 2019
### TMEM57 -> MACO1
Songyang.25<-c('ABAT', 'BCAR3', 'CTSF', 'DEAF1', 'ENC1', 'ETV5', 'FAM117A', 'FZD2', 'GALNT12', 'GALNT3', 'GJB3',
               'KDM6A', 'KYNU', 'PCNA', 'PFKP', 'PLEK2', 'RASGRP2', 'SERPIND1', 'SGSH',
               'TLE1', 'TMEM38B', 'MACO1', 'TRIM45', 'USP47', 'VWA1')

## Zhang et al. 2019
## ERO1L -> ERO1A
Zhang.4<-c('MYO1E', 'ERO1A', 'C1QTNF6', 'FAM83A')

## Yue et al. 2019
Yue.3<-c('ADAM12', 'BTK', 'ERG')

## Xie. 2019
Xie.6<-c('RRAGB', 'RSPH9', 'RPS6KL1', 'RXFP1', 'RRM2', 'RTL1')

## Jiawei et al. 2020
Jiawei.6<-c('VIPR1', 'FCN3', 'CA4', 'CYP4B1', 'CRTAC1', 'NEDD9')

## Liu et al. 2020
Liu.14<-c('C1QTNF6', 'ERO1A', 'MELTF', 'ITGB1-DT', 'RGS20', 'FETUB', 'NTSR1', 'LINC02178', 'AC034223.2', 'LINC01312',
          'AL353746.1', 'AC034223.1', 'DRAXINP1', 'LINC02310')

## Zengin et al. 2020
Zengin.12<-c('BCHE', 'CCNA1', 'CYP24A1', 'DEPTOR', 'MASP2', 'MGLL', 'MYO1A', 'PODXL2', 'RAPGEF3', 'SGK2', 'TNNI2', 'ZBTB16')

## Wu et al. 2020
Wu.21<-c('USP7', 'SPHK1', 'SMAD6', 'RIPK2', 'RAC1', 'PTCH1', 'PMAIP1', 'PLAUR', 'MOV10', 'MMP12', 'MIF', 'ITPR1', 'IL6ST', 
         'IL32', 'HMOX1', 'ELF4', 'C7', 'C5AR1', 'BIRC5', 'ARF6', 'AQP3')

## Liu et al. 2021
Liu.4<-c('CENPH', 'MYLIP', 'PITX3', 'TRAF3IP3')

## Al-Dherasi et al. 2021
Al_Dherasi.7<-c('UCN2', 'RIMS2', 'CAVIN2', 'GRIA1', 'PKHD1L1', 'PGM5' , 'CLIC6')

## Wu et al. 2022
Wu.4<-c('HLF', 'CHRDL1', 'SELENBP1', 'TMEM163')

## Sun et al. 2022
Sun.14<-c('CPS1', 'CTPS2', 'DARS2', 'IGFBP3', 'MCM5', 'MCM7', 'NME4', 'NT5E',
          'PLK1', 'POLR3G', 'PTTG1', 'SERPINB5', 'TXNRD1', 'TYMS')

## Dessie et al. 2022
Dessie.9<-c('C1QTNF6', 'CDC25C', 'E2F7', 'IL11', 'CLEC12B', 'CYP17A1', 'FAM72D', 'GNG7', 'RTN1')

## Zhu et al. 2023
Zhu.10<-c('ACAN', 'ADAMTS15', 'ADAMTS8', 'BCAN', 'COL4A3', 'ITGA8', 'ITGB4', 'LAD1', 'TENM3', 'TIMP1')

## Xia et al. 2023
Xia.5<-c('TCN1', 'CENPF', 'MAOB', 'CRTAC1','PLEK2')



### previous marker added (EGFR-specific markers)
Li.3<-c('BTLA', 'BUB1B', 'CENPE')
Zhang.6<-c("B3GNT3", 'CDH3', 'CST1', 'KLB', 'KRT15', 'ZBTB16')

previous_markers<-list(CIN25, Zhu.15, Park.2012, CCP.markers, Lu.16, Okayama.4, Gentles.9, ESLA7, Li.8,
                       Songyang.25, Zhang.4, Yue.3, Xie.6, Jiawei.6, Zengin.12, Wu.21, Liu.4, Al_Dherasi.7, Wu.4, Sun.14,
                       Dessie.9, Zhu.10, Xia.5, Li.3, Zhang.6)
names(previous_markers)<-c('Carter_2006', 'Zhu_2010', 'Park_2012', 'Wistuba_2013', 'Lu_2013',
                           'Okayama_2014', 'Gentles_2015', 'Krzystanek_2016', 'Li_2018', 'Songyang_2019', 
                           'Zhang_2019', 'Yue_2019', 'Xie_2019', 'Jiawei_2020', 'Zengin_2020', 'Wu_2020',
                           'Liu_2021', 'Al-Dherasi_2021', 'Wu_2022',
                           'Sun_2022', 'Dessie_2022', 'Zhu_2023', 'Xia_2023', 'Li_2021', 'Zhang_2021')

