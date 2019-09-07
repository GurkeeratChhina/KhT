field = 2

name1 = "PretzelTangles/(2,-3)"
Tangle1 = Tangle.PretzelTangle(2, -3)
BNcx1 = Tangle1.toReduced_BNComplex(1000, field)
multicurve1 = BNcx1.to_multicurve()
multicurve1.drawpng(name1+"_BNr","hdelta",Tangle1.slices)