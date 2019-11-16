field = 2

html_content=""

for i in range(1,6):
    for j in range(1,6):
        name = "("+str(2*i)+","+str(-2*j-1)+")-PretzelTangle"
        Tangle = Tangle.PretzelTangle(2*i, -2*j-1)
        BNcx = Tangle.toReduced_BNComplex(1000, field)
        BNcx = BNcx.cone(1)
        BNcx.clean_up(1000)
        multicurve = BNcx.to_multicurve()
        multicurve.save(name)
        html_content += multicurve.html(name,"_BNr","hdelta",Tangle)
          
header="""---
title: Pretzel Tangles
layout: default
filename: Pretzel_Tangle
---
<link rel="stylesheet" href="css/main.css">
<script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
<script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>

<script>
function on(file,pagenumber) {
    document.getElementById("overlay_content").src=String(file)+".pdf";
    document.getElementById("overlay_title").innerHTML=String(file);
    document.getElementById("overlay").style.display = "block";
}

function off() {
    document.getElementById("overlay").style.display = "none";
} 
</script>

<div id="overlay" onclick="off()" >
<div id="text">
<div id="overlay_title_frame"><h4 id="overlay_title"></h4></div>
<iframe id="overlay_content" width="500px" height="500px">
</iframe></div>
</div> 

<h4>Reduced arc-invariant over \(\mathbb{F}_2\)</h4>
"""

with open("examples/"+filename+".html", "w") as text_file:
    print(header+html_content, file=text_file)
