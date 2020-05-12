field = 2

html_content=""

name = "TorusKnot"
#Tangle = Tangle("cap1.cap2.cap3.neg2.neg1.neg0.neg2.neg1.neg0.neg2.neg1.neg0.neg2.neg1.neg0.neg2.neg1.neg0.cup3.cup2.neg0.neg0")
Tangle = Tangle("cap1.cap2.cap3.neg2.neg1.neg0.neg2.neg1.neg0.neg2.neg1.neg0.neg2.neg1.neg0.neg2.neg1.neg0.cup3.cup2")
BNcx = Tangle.toReduced_BNComplex(10000, field)
#BNcx = BNcx.cone(1)
BNcx.clean_up(10000)
multicurve = BNcx.to_multicurve()
multicurve.save(name)
multicurve.html(name,"_BNr_F2","hdelta",Tangle)

