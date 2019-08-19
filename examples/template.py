name = "template" 
tangle = "cup2.cup3.neg1.neg2.neg0.neg1.pos3.pos4.pos2.pos3.neg1.neg0.neg2.neg1.cap3.cap2.cap1" # tangle data
drawtangle(tangle,name,"slices",1) # draw tangle

# cx = BNbracket(tangle,0,0,1) # compute Bar-Natan's bracket
# BNr = cx.ToBNAlgebra(7) # convert the Bar-Natan's bracket into a complex over BNAlgebra

# cx.save(name) # save the complex cx for later
# cx = importCobcx(name) # import saved complex

# complex1 = BN_complex.ToCob() # conversion back to Cob does not work yet, since Z is not implemented over BNAlgebra and Cob is only implemented over Z.

# BNr.eliminateAll() # cancel all identity components of the differential
# BNr.clean_up() # try to find the immersed curve invariant BNr through a sequence of random isotopies
# multicurve = BNr.to_multicurve() # if the previous step was successful, separate the components
# multicurve.draw(name+"_BNr","hdelta",tangle) # create output pdf-file

# Khr=BNr.cone(1) # compute the complex Khr as a cone [BNr---H¹--->BNr]
# Khr.clean_up() # try to find the immersed curve invariant Khr through a sequence of random isotopies
# multicurve_KhT=Khr.to_multicurve() # if the previous step was successful, separate the components
# multicurve_KhT.draw(name+"_Khr","hdelta",tangle) # create output pdf-file

# BNr.save(name) # save this complex for later
# BNr = importBNcx(name) # import saved complex
