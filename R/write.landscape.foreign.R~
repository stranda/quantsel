"landscape.write.foreign" <-
function(rland, numi=24, fn = "foreign", fmt="GDA")
{
  if  (is.landscape(rland))
    {
      if (fmt %in% c("GDA","Gda","gda"))
        {
          .Call("writeGDA",fn,rland,numi,PACKAGE = "kernelPop2")
        }
      if (fmt %in% c("Arlequin","arlequin","ArlequinDip","Arlequindip","arlequindip"))
        {
          .Call("writeArlequinDip",fn,rland,numi,PACKAGE = "kernelPop2")
        }
      if (fmt %in% c("ArlequinHap","Arlequinhap","arlequinhap"))
        {
          print("arlequin haploid export engine not currently working")
#          .Call("writeArlequinHap",fn,rland,numi,PACKAGE = "kernelPop2")
        }
      if (fmt %in% c("BIOSYS-1","BIOSYS","Biosys","biosys","biosys-1"))
        {
          .Call("writeBIOSYS",fn,rland,numi,PACKAGE = "kernelPop2")
        }
      if (fmt %in% c("GenPop","Genpop","genpop"))
        {
          print("a better alternative is landscape.genepop.output()")
          .Call("writeGenPop",fn,rland,numi,PACKAGE = "kernelPop2")
        }
      if (fmt %in% c("R","r"))
        {
          .Call("writeR",fn,rland,numi,PACKAGE = "kernelPop2")
        }
      if (fmt %in% c("Migrate","MigrateDiploid","migrate","migratediploid"))
        {
          .Call("writeMigrateDip",fn,rland,numi,PACKAGE = "kernelPop2")
        }
      if (fmt %in% c("ReRat","rerat"))
        {
          .Call("writeGenPop",fn,rland,numi,PACKAGE = "kernelPop2")
        }
    }
  else
    {
      print ("rland not in landscape format")
    }
}

