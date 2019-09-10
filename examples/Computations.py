field = 2

name20 = "../Output/LiamsTangle_2"
Tangle20 = Tangle.LiamsTangle(2, [0, 0])
BNcx20 = Tangle20.toReduced_BNComplex(1000, field)
multicurve20 = BNcx20.to_multicurve()
multicurve20.draw(name20+"_BNr","hdelta",Tangle20.slices)

name21 = "../Output/LiamsTangle_3"
Tangle21 = Tangle.LiamsTangle(3, [0, 0, 0])
BNcx21 = Tangle21.toReduced_BNComplex(1000, field)
multicurve21 = BNcx21.to_multicurve()
multicurve21.draw(name21+"_BNr","hdelta",Tangle21.slices)

name22 = "../Output/LiamsTangle_1_[1]"
Tangle22 = Tangle.LiamsTangle(1, [1])
BNcx22 = Tangle22.toReduced_BNComplex(1000, field)
multicurve22 = BNcx22.to_multicurve()
multicurve22.draw(name22+"_BNr","hdelta",Tangle22.slices)

name23 = "../Output/LiamsTangle_1_[2]"
Tangle23 = Tangle.LiamsTangle(1, [2])
BNcx23 = Tangle23.toReduced_BNComplex(1000, field)
multicurve23 = BNcx23.to_multicurve()
multicurve23.draw(name23+"_BNr","hdelta",Tangle23.slices)

name24 = "../Output/LiamsTangle_1_[3]"
Tangle24 = Tangle.LiamsTangle(1, [3])
BNcx24 = Tangle24.toReduced_BNComplex(1000, field)
multicurve24 = BNcx24.to_multicurve()
multicurve24.draw(name24+"_BNr","hdelta",Tangle24.slices)

name25 = "../Output/LiamsTangle_2_[1,1]"
Tangle25 = Tangle.LiamsTangle(2, [1, 1])
BNcx25 = Tangle25.toReduced_BNComplex(1000, field)
multicurve25 = BNcx25.to_multicurve()
multicurve25.draw(name25+"_BNr","hdelta",Tangle25.slices)

name26 = "../Output/LiamsTangle_2_[1,-1]"
Tangle26 = Tangle.LiamsTangle(2, [1, -1])
BNcx26 = Tangle26.toReduced_BNComplex(1000, field)
multicurve26 = BNcx26.to_multicurve()
multicurve26.draw(name26+"_BNr","hdelta",Tangle26.slices)

name27 = "../Output/LiamsTangle_2_[1,-2]"
Tangle27 = Tangle.LiamsTangle(2, [1, -2])
BNcx27 = Tangle27.toReduced_BNComplex(1000, field)
multicurve27 = BNcx27.to_multicurve()
multicurve27.draw(name27+"_BNr","hdelta",Tangle27.slices)

name28 = "../Output/LiamsTangle_2_[1,2]"
Tangle28 = Tangle.LiamsTangle(2, [1, 2])
BNcx28 = Tangle28.toReduced_BNComplex(1000, field)
multicurve28 = BNcx28.to_multicurve()
multicurve28.draw(name28+"_BNr","hdelta",Tangle28.slices)

name29 = "../Output/LiamsTangle_2_[2,-2]"
Tangle29 = Tangle.LiamsTangle(2, [2, -2])
BNcx29 = Tangle29.toReduced_BNComplex(1000, field)
multicurve29 = BNcx29.to_multicurve()
multicurve29.draw(name29+"_BNr","hdelta",Tangle29.slices)

name30 = "../Output/LiamsTangle_2_[2,2]"
Tangle30 = Tangle.LiamsTangle(2, [2, 2])
BNcx30 = Tangle30.toReduced_BNComplex(1000, field)
multicurve30 = BNcx30.to_multicurve()
multicurve30.draw(name30+"_BNr","hdelta",Tangle30.slices)

