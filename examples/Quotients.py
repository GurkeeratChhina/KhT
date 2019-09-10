field = 5

name1 = "Quotients_of_(2,-3)_PretzelTangle/NoTwists"
Tangle1 = Tangle.quotient_of_2_m3_pretzel_tangle(0)
BNcx1 = Tangle1.toReduced_BNComplex(1000, field)
multicurve1 = BNcx1.to_multicurve()
multicurve1.drawpng(name1 + "_BNr", "hdelta", Tangle1.slices)

name2 = "Quotients_of_(2,-3)_PretzelTangle/+1Twist"
Tangle2 = Tangle.quotient_of_2_m3_pretzel_tangle(1)
BNcx2 = Tangle2.toReduced_BNComplex(1000, field)
multicurve2 = BNcx2.to_multicurve()
multicurve2.drawpng(name2 + "_BNr", "hdelta", Tangle2.slices)

name3 = "Quotients_of_(2,-3)_PretzelTangle/+2Twists"
Tangle3 = Tangle.quotient_of_2_m3_pretzel_tangle(2)
BNcx3 = Tangle3.toReduced_BNComplex(1000, field)
multicurve3 = BNcx3.to_multicurve()
multicurve3.drawpng(name3 + "_BNr", "hdelta", Tangle3.slices)

name4 = "Quotients_of_(2,-3)_PretzelTangle/+3Twists"
Tangle4 = Tangle.quotient_of_2_m3_pretzel_tangle(3)
BNcx4 = Tangle4.toReduced_BNComplex(1000, field)
multicurve4 = BNcx4.to_multicurve()
multicurve4.drawpng(name4 + "_BNr", "hdelta", Tangle4.slices)

name5 = "Quotients_of_(2,-3)_PretzelTangle/+4Twists"
Tangle5 = Tangle.quotient_of_2_m3_pretzel_tangle(4)
BNcx5 = Tangle5.toReduced_BNComplex(1000, field)
multicurve5 = BNcx5.to_multicurve()
multicurve5.drawpng(name5 + "_BNr", "hdelta", Tangle5.slices)

name6 = "Quotients_of_(2,-3)_PretzelTangle/-1Twist"
Tangle6 = Tangle.quotient_of_2_m3_pretzel_tangle(-1)
BNcx6 = Tangle6.toReduced_BNComplex(1000, field)
multicurve6 = BNcx6.to_multicurve()
multicurve6.drawpng(name6 + "_BNr", "hdelta", Tangle6.slices)

name7 = "Quotients_of_(2,-3)_PretzelTangle/-2Twists"
Tangle7 = Tangle.quotient_of_2_m3_pretzel_tangle(-2)
BNcx7 = Tangle7.toReduced_BNComplex(1000, field)
multicurve7 = BNcx7.to_multicurve()
multicurve7.drawpng(name7 + "_BNr", "hdelta", Tangle7.slices)