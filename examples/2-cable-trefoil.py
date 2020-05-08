name = "2-cable-trefoil" 
tangle = "cap1.cap2.cap3.neg1.neg2.neg0.neg1.pos3.pos2.pos4.pos3.neg1.neg0.neg2.neg1.cup3.cup2"
Tangle = Tangle(tangle)
html_content=""

#cx.save(name) # save this complex for later
#cx = importCobcx(name) # import saved complex
#complex1 = BN_complex.ToCob() # conversion back to Cob does not work yet, since Z is not implemented over BNAlgebra and Cob is only implemented over Z.

#drawtangle(tangle,name+"test",start=1,title=[name+"test",""])

#Khr=BNr.cone(1) # compute the complex Khr as a cone [BNr---H¹--->BNr]
#Khr.clean_up() # try to find the immersed curve invariant Khr through a sequence of random isotopies
#multicurve_KhT=Khr.to_multicurve() # if the previous step was successful, separate the components
#multicurve_KhT.draw(name+"_Khr","hdelta",tangle,thumbnails=True) # create output pdf-file

#BNr.save(name) # save this complex for later
#BNr = importBNcx(name) # import saved complex

html_content += "<h4>Arc invariant over \(\mathbb{F}_2\)</h4>"
cx = BNbracket(tangle,0,0,1) # compute Bar-Natan's bracket
BNr = cx.ToBNAlgebra(2) # convert the Bar-Natan's bracket into a complex over BNAlgebra
BNr.eliminateAll() # cancel all identity components of the differential
BNr.clean_up() # try to find the immersed curve invariant BNr through a sequence of random isotopies
#BNr.draw(name)

multicurve = BNr.to_multicurve()
multicurve.save(name)
html_content += multicurve.html(name,"_BNr7","hdelta",Tangle)

html_content += "<h4>Figure-8 invariant over \(\mathbb{F}_7\)</h4>"

Khr = BNr.cone(1)
Khr.clean_up()
multicurveKh = Khr.to_multicurve()
html_content += multicurveKh.html(name,"_Khr7","hdelta",Tangle)

html_content += "<h4>Arc invariant over \(\mathbb{F}_2\)</h4>"
cx = BNbracket(tangle,0,0,1) # compute Bar-Natan's bracket
BNr = cx.ToBNAlgebra() # convert the Bar-Natan's bracket into a complex over BNAlgebra
BNr.eliminateAll() # cancel all identity components of the differential
BNr.clean_up() # try to find the immersed curve invariant BNr through a sequence of random isotopies
#BNr.draw(name)

multicurve = BNr.to_multicurve()
multicurve.save(name)
html_content += multicurve.html(name,"_BNr2","hdelta",Tangle)

html_content += "<h4>Figure-8 invariant over \(\mathbb{F}_2\)</h4>"

Khr = BNr.cone(1)
Khr.clean_up()
multicurveKh = Khr.to_multicurve()
html_content += multicurveKh.html(name,"_Khr2","hdelta",Tangle)
        
          
header="""---
title: 2-cable trefoil
layout: default
filename: 2-cable-trefoil
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
"""

with open("examples/"+filename+".html", "w") as text_file:
    print(header+html_content, file=text_file)
