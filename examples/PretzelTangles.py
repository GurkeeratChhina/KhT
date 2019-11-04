field = 2

content=""

for i in range(1,6):
    for j in range(1,6):
        name = "PretzelTangles-("+str(2*i)+","+str(-2*j-1)+")"
        Tangle = Tangle.PretzelTangle(2*i, -2*j-1)
        BNcx = Tangle.toReduced_BNComplex(1000, field)
        multicurve = BNcx.to_multicurve()
        multicurve.draw(name+"_BNr_thumbs","hdelta",[name+"_BNr",""],Tangle.slices,thumbnails=True)
        multicurve.draw(name+"_BNr","qhdelta",[name+"_BNr",""],Tangle.slices)
        
        # convert thumbnails file into pngs
        run("cd examples && pdftoppm '"+name+"_BNr_thumbs.pdf' '"+name+"_BNr_thumbs' -png -rx 300 -ry 300", shell=True)
        # split detailed pdf into single pages
        run("cd examples && pdftk '"+name+"_BNr.pdf' burst output '"+name+"_BNr'_%02d.pdf", shell=True)  
        
        content+="""
        <tr>"""
        for counter in range(len(multicurve.comps)+1):
            counterstring=str(counter+1)
            if len(counterstring)==1:
                counterstring="0"+counterstring
            content+="""<td><a href="{0}_BNr_{2}.pdf"><img src="{0}_BNr_thumbs-{1}.png" style="height:200px;"></a> </td>""".format(name,counter+1,counterstring)
        
        content+="""
        </tr>
        """
header="""---
title: Pretzel Tangles
layout: default
filename: Pretzel_Tangle
---

<script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
<script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>

<table>
    <thead>
        <tr>
            <th>Tangle</th>
	    <th colspan = "4" >reduced arc-invariant over \(\mathbb{F}_2\)</th>
        </tr>
    </thead>
    <tbody>
"""
footer="""
    </tbody>
</table>
"""

with open("examples/pretzels.html", "w") as text_file:
    print(header+content+footer, file=text_file)
