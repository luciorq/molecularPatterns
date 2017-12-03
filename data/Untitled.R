
#test  = loadSam("data/SPO-108.VSVgfp.v1.sam",15,30)

#test distribution and density
#utils::head(test)
#p=plotSizeDistribution(samObject = test,norm = 30000)
#p=plotSizeDistribution(samObject = test,norm = 30000,sizeStart = 18,sizeEnd = 25)
#calcDensityPerBase(test,"VSVgfp_Olmo",ymin=-10000,ymax=10000)

#ACF
#samObject =test
#calcACF(test, ref)
#ref='KU721836_synthetic'
#regIni = 1
#regEnd=1000
#sizeStart=24
#sizeEnd=24
#
##weblogo
#samObject =test
#ref='KU721836_synthetic'
#sizeStart =26
#sizeEnd   =30
#createWebLogo(samObject, ref)
#
#
#
##phasing / offset test
#test  = loadSam("data/SPO-108_10.synthetic_VSV.sam",24,30)
#test  = loadSamLikeBed("data/teste_h100.sam",15,30)
#samObject =test
#nrow(test)
#ref='KU721836_synthetic'
#extend_upstream=20
#extend_downstream=20
#readSize=0
#bin=1
#min=1
#max=12015
#utils::head(samObject)
#
