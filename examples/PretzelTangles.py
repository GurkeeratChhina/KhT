field = 2

html_content=""

for i in range(1,3):
    for j in range(1,3):
        name = "("+str(2*i)+","+str(-2*j-1)+")-PretzelTangle"
        Tangle = Tangle.PretzelTangle(2*i, -2*j-1)
        #if i<=j:
        #    Tangle.twist([-2*i])
        #elif i==j+1:# arc is a rational component of slope 1/0
        #    Tangle.twist([2*j+2])
        #else:# i>j+1 # add 2*(i-j)-1 twists on the right
        #    Tangle.twist([2*j+1])
        #    #Tangle.twist([2*j+1,1-2*(i-j)])
        BNcx = Tangle.toReduced_BNComplex(10000, field)
        #BNcx = BNcx.cone(1)
        BNcx.clean_up(10000)
        multicurve = BNcx.to_multicurve()
        multicurve.save(name)
        html_content+="(i,j)=("+str(i)+","+str(j)+"):"
        html_content += multicurve.html(name,"_BNr","hdelta",Tangle)

html_content += "<h4>Figure-8 invariant over \(\mathbb{F}_2\)</h4>"

for i in range(1,3):
    for j in range(1,3):
        name = "("+str(2*i)+","+str(-2*j-1)+")-PretzelTangle"
        Tangle = Tangle.PretzelTangle(2*i, -2*j-1)
        #if i<=j:
        #    Tangle.twist([-2*i])
        #elif i==j+1:
        #    Tangle.twist([2*j+2])
        #else:
        #    Tangle.twist([2*j+1])
        #    #Tangle.twist([2*j+1,1-2*(i-j)])
        BNcx = Tangle.toReduced_BNComplex(10000, field)
        BNcx = BNcx.cone(1)
        BNcx.clean_up(1000)
        multicurve = BNcx.to_multicurve()
        multicurve.save(name)
        html_content+=str([i,j])
        html_content += multicurve.html(name,"_Khr","hdelta",Tangle)
          
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

<h4>Arc invariant over \(\mathbb{F}_2\)</h4>
"""

with open("examples/"+filename+".html", "w") as text_file:
    print(header+html_content, file=text_file)
