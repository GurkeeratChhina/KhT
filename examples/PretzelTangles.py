field = 2

for i in range(1,5):
    for j in range(1,5):
        name = "PretzelTangles-("+str(2*i)+","+str(-2*j-1)+")"
        Tangle = Tangle.PretzelTangle(2*i, -2*j-1)
        BNcx = Tangle.toReduced_BNComplex(1000, field)
        multicurve = BNcx.to_multicurve()
        multicurve.draw(name+"_BNr_thumbs","hdelta",[name+"_BNr",""],Tangle.slices,thumbnails=True)
        multicurve.draw(name+"_BNr","qhdelta",[name+"_BNr",""],Tangle.slices)
