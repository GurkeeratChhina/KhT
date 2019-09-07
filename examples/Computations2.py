field = 2

name1 = "../Output/trefoil_linked_with_line"
Tangle1 = Tangle("cup0.pos1.neg2.pos1.neg0.neg1.pos2.cap3.cap1")
BNcx1 = Tangle1.toReduced_BNComplex(1000, field)
multicurve1 = BNcx1.to_multicurve()
multicurve1.draw(name1+"_BNr","hdelta",Tangle1.slices)

name2 = "../Output/pseudo_trefoil_linked_with_line"
Tangle2 = Tangle("cup0.neg1.pos0.pos2.neg1.pos2.pos0.cap3.cap1")
BNcx2 = Tangle2.toReduced_BNComplex(1000, field)
multicurve2 = BNcx2.to_multicurve()
multicurve2.draw(name2+"_BNr","hdelta",Tangle2.slices)

name3 = "../Output/twist_linked_with_cap"
Tangle3 = Tangle("cup1.pos2.pos0.pos0.neg1.neg1.pos2.neg1.pos0.cap3.cap1")
BNcx3 = Tangle3.toReduced_BNComplex(1000, field)
multicurve3 = BNcx3.to_multicurve()
multicurve3.draw(name3+"_BNr","hdelta",Tangle3.slices)

name4 = "../Output/pretzel_3_minus_3_looped_bottom"
Tangle4 = Tangle("cup1.pos0.neg1.neg2.neg2.neg2.pos0.pos0.pos0.cap3.cap1")
BNcx4 = Tangle4.toReduced_BNComplex(1000, field)
multicurve4 = BNcx4.to_multicurve()
multicurve4.draw(name4+"_BNr","hdelta",Tangle4.slices)

name5 = "../Output/cable-trefoil"
Tangle5 = Tangle("cup1.neg0.pos1.neg0.cap1", None,1,1).Cable()
BNcx5 = Tangle5.toReduced_BNComplex(1000, field)
multicurve5 = BNcx5.to_multicurve()
multicurve5.draw(name5+"_BNr","hdelta",Tangle5.slices)

name6 = "../Output/cable-fig8"
Tangle6 = Tangle("cup1.neg0.neg0.pos1.neg0.cap1", None,1,1).Cable()
BNcx6 = Tangle6.toReduced_BNComplex(1000, field)
multicurve6 = BNcx6.to_multicurve()
multicurve6.draw(name6+"_BNr","hdelta",Tangle6.slices)

name7 = "../Output/pretzel_2_minus_2_looped_bottom"
Tangle7 = Tangle("cup1.pos0.neg1.neg2.neg2.pos0.pos0.cap3.cap1")
BNcx7 = Tangle7.toReduced_BNComplex(1000, field)
multicurve7 = BNcx7.to_multicurve()
multicurve7.draw(name7+"_BNr","hdelta",Tangle7.slices)

name8 = "../Output/LiamsTangle"
Tangle8 = Tangle.LiamsTangle(1, [0])
BNcx8 = Tangle8.toReduced_BNComplex(1000, field)
multicurve8 = BNcx8.to_multicurve()
multicurve8.draw(name8+"_BNr","hdelta",Tangle8.slices)

name9 = "../Output/two_twist_hitch"
Tangle9 = Tangle("cup1.cup1.pos4.neg0.neg3.pos1.cap2.pos2.pos2.neg0.neg0.cap3.cap1")
BNcx9 = Tangle9.toReduced_BNComplex(1000, field)
multicurve9 = BNcx9.to_multicurve()
multicurve9.draw(name9+"_BNr","hdelta",Tangle9.slices)

name10 = "../Output/quotient_of_2_m3_pretzel_tangle"
Tangle10 = Tangle.quotient_of_2_m3_pretzel_tangle(0)
BNcx10 = Tangle10.toReduced_BNComplex(1000, field)
multicurve10 = BNcx10.to_multicurve()
multicurve10.draw(name10 + "_BNr", "hdelta", Tangle10.slices)

name11 = "../Output/cable-5_1"
Tangle11 = Tangle("cup1.pos0.neg1.neg1.neg1.pos0.cap1", None,1,1).Cable()
BNcx11 = Tangle11.toReduced_BNComplex(1000, field)
multicurve11 = BNcx11.to_multicurve()
multicurve11.draw(name11+"_BNr","hdelta",Tangle11.slices)

name12 = "../Output/cable-5_2"
Tangle12 = Tangle("cup1.neg0.pos1.pos1.neg0.neg0.cap1", None,1,1).Cable()
BNcx12 = Tangle12.toReduced_BNComplex(1000, field)
multicurve12 = BNcx12.to_multicurve()
multicurve12.draw(name12+"_BNr","hdelta",Tangle12.slices)

name13 = "../Output/cable-6_1"
Tangle13 = Tangle("cup1.neg0.pos1.pos1.pos1.neg0.neg0.cap1", None,1,1).Cable()
BNcx13 = Tangle13.toReduced_BNComplex(1000, field)
multicurve13 = BNcx13.to_multicurve()
multicurve13.draw(name13+"_BNr","hdelta",Tangle13.slices)

name40 = "../Output/quotient_of_2_m3_pretzel_tangle_minus_1"
Tangle40 = Tangle.quotient_of_2_m3_pretzel_tangle(-1)
BNcx40 = Tangle40.toReduced_BNComplex(1000, field)
multicurve40 = BNcx40.to_multicurve()
multicurve40.draw(name40 + "_BNr", "hdelta", Tangle40.slices)

name41 = "../Output/quotient_of_2_m3_pretzel_tangle_minus_2"
Tangle41 = Tangle.quotient_of_2_m3_pretzel_tangle(-2)
BNcx41 = Tangle41.toReduced_BNComplex(1000, field)
multicurve41 = BNcx41.to_multicurve()
multicurve41.draw(name41 + "_BNr", "hdelta", Tangle41.slices)
