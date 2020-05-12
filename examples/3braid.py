field = 2
name = "3-braid"
Tangle = Tangle("cap1.cap3."+("neg2."*21)+("pos0."*20)+"cup1")
#Tangle = Tangle.PretzelTangle(2*i, -20)
#r1.l3.x2.x2.x2.x2.x2.x2.x2.x2.x2.x2.x2.y0.y0.y0.y0.y0.y0.y0.y0.y0.y0.u1
BNcx = Tangle.toReduced_BNComplex(10000, field)

