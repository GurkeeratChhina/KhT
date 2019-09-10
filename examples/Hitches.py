field = 5

name1 = "TwoTwistHitch/NoCrossings"
Tangle1 = Tangle.two_twist_hitch(0)
BNcx1 = Tangle1.toReduced_BNComplex(1000, field)
multicurve1 = BNcx1.to_multicurve()
multicurve1.drawpng(name1+"_BNr","hdelta",Tangle1.slices)

name2 = "TwoTwistHitch/+1Crossing"
Tangle2 = Tangle.two_twist_hitch(1)
BNcx2 = Tangle2.toReduced_BNComplex(1000, field)
multicurve2 = BNcx2.to_multicurve()
multicurve2.drawpng(name2+"_BNr","hdelta",Tangle2.slices)

name3 = "TwoTwistHitch/+2Crossing"
Tangle3 = Tangle.two_twist_hitch(2)
BNcx3 = Tangle3.toReduced_BNComplex(1000, field)
multicurve3 = BNcx3.to_multicurve()
multicurve3.drawpng(name3+"_BNr","hdelta",Tangle3.slices)

name4 = "TwoTwistHitch/+3Crossing"
Tangle4 = Tangle.two_twist_hitch(3)
BNcx4 = Tangle4.toReduced_BNComplex(1000, field)
multicurve4 = BNcx4.to_multicurve()
multicurve4.drawpng(name4+"_BNr","hdelta",Tangle4.slices)

name5 = "TwoTwistHitch/+4Crossing"
Tangle5 = Tangle.two_twist_hitch(4)
BNcx5 = Tangle5.toReduced_BNComplex(1000, field)
multicurve5 = BNcx5.to_multicurve()
multicurve5.drawpng(name5+"_BNr","hdelta",Tangle5.slices)

name6 = "TwoTwistHitch/-1Crossing"
Tangle6 = Tangle.two_twist_hitch(-1)
BNcx6 = Tangle6.toReduced_BNComplex(1000, field)
multicurve6 = BNcx6.to_multicurve()
multicurve6.drawpng(name6+"_BNr","hdelta",Tangle6.slices)

name7 = "TwoTwistHitch/-2Crossing"
Tangle7 = Tangle.two_twist_hitch(-2)
BNcx7 = Tangle7.toReduced_BNComplex(1000, field)
multicurve7 = BNcx7.to_multicurve()
multicurve7.drawpng(name7+"_BNr","hdelta",Tangle7.slices)