# -*- coding: utf-8 -*-
# COPYRIGHT 2019 Gurkeerat Chhina, Claudius Zibrowius
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import cairo
from subprocess import run

# The drawing code is not used often - now we use the complex drawing code. 
    
def SanitizeTeX(string):
    return string.replace("_", "\_")

def draw_tangle_ends(posx,posy,clt,h,ctx):
    ctx.set_font_size(0.40)
    ctx.select_font_face("Courier",cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_BOLD)
    
    #dots with labels
    for i in range(0,clt.top):
        ctx.move_to(i+posx,posy)
        ctx.line_to(i+posx,posy)
        s = '{}'.format(i)
        xbearing, ybearing, width, height, dx, dy = ctx.text_extents(s)
        ctx.move_to(i - width/2+posx,posy-0.15)
        ctx.show_text(s)
    
    for i in range(0,clt.bot):
        ctx.move_to(i+posx,h+posy)
        ctx.line_to(i+posx,h+posy)
        s = '{}'.format(i+clt.top)
        xbearing, ybearing, width, height, dx, dy = ctx.text_extents(s)
        ctx.move_to(i - width/2+posx,posy+h+0.35)
        ctx.show_text(s)
    
    ctx.set_source_rgb(0,0,0)
    ctx.set_line_width(0.1)
    ctx.stroke()
    
def draw_arcs(posx,posy,clt,h,ctx,rgb):
    
    #arcs connecting dots
    for i in clt.arcs_top():
        ctx.move_to(posx+i[0],posy+0)
        ctx.curve_to(posx+i[0],posy+abs(i[1]-i[0])/2,posx+i[1],posy+abs(i[1]-i[0])/2,posx+i[1],posy+0)
    
    for i in clt.arcs_bot():
        ctx.move_to(posx+i[0]-clt.top,posy+h)
        ctx.curve_to(posx+i[0]-clt.top,posy+h-abs(i[1]-i[0])/2,posx+i[1]-clt.top,posy+h-abs(i[1]-i[0])/2,posx+i[1]-clt.top,posy+h)
        
    for i in clt.arcs_mix():
        ctx.move_to(posx+i[0],posy+0)
        ctx.curve_to(posx+i[0],posy+abs(i[1]-clt.top-i[0])/2,posx+i[1]-clt.top,posy+h-abs(i[1]-clt.top-i[0])/2,posx+i[1]-clt.top,posy+h)
    
    ctx.set_source_rgb(*rgb)
    ctx.set_line_width(0.06)
    ctx.stroke()
    
def draw_dot_on_arc(arc,clt,h,ctx,deco_index):
    ctx.set_line_width(0.3)
    
    if arc[0]<clt.top:
        ctx.move_to(2+arc[0],(h+1)*deco_index+0.1)
        ctx.line_to(2+arc[0],(h+1)*deco_index+0.1)
    else:
        ctx.move_to(2+arc[0]-clt.top,h+(h+1)*deco_index-0.1)
        ctx.line_to(2+arc[0]-clt.top,h+(h+1)*deco_index-0.1)
    
    # THE CODE IS NOT OPERATIONAL YET. SO FAR, WE DRAW COLOURED DOTS ON TANGLE ENDS, BUT IDEALLY DRAW DOTS ON ARCS. 
    #arcs connecting dots
    #if arc not in clt.arcs():
    #    raise Exception('Proper dot-decoration of cobordisms failed.')
    # 
    # if arc in clt.arcs_top():
    #    ctx.move_to(posx+i[0],posy+0)
    #    ctx.curve_to(posx+i[0],posy+abs(i[1]-i[0])/2,posx+i[1],posy+abs(i[1]-i[0])/2,posx+i[1],posy+0)
    #
    #if arc in clt.arcs_bot():
    #    ctx.move_to(posx+i[0]-clt.top,posy+h)
    #    ctx.curve_to(posx+i[0]-clt.top,posy+h-abs(i[1]-i[0])/2,posx+i[1]-clt.top,posy+h-abs(i[1]-i[0])/2,posx+i[1]-clt.top,posy+h)
    #    
    #if arc in clt.arcs_mix():
    #    ctx.move_to(posx+i[0],posy+0)
    #    ctx.curve_to(posx+i[0],posy+abs(i[1]-clt.top-i[0])/2,posx+i[1]-clt.top,posy+h-abs(i[1]-clt.top-i[0])/2,posx+i[1]-clt.top,posy+h)
    #
    ctx.set_source_rgb(255,0,0)
    ctx.stroke()

def drawtangle(string,filename,style="plain",start=1,subtitle=""):
    """draw the tangle specified by 'string', which is a concatenation of words <type>+<index>, separated by '.' read from right to left, for each elementary tangle slice, read from top to bottom, where:
    <type> is equal to:
        'pos': positive crossing
        'neg': negative crossing
        'cup': cup
        'cap': cap
    <index> is the index at which the crossing, cap or cup sits. 
    The optional parameter 'style' is either 'plain' (simple tangle) or 'slices' (shows slices with labels).
    The optional parameter 'start' is an integer which specifies the number of tangle ends at the top.
    The optional parameter 'subtitle' is a string which is added to the top.
    """
    content="""\\documentclass{article}
\\usepackage[crop=off]{auto-pst-pdf}
\\usepackage{pst-node,rotating}
\\renewcommand{\\familydefault}{\\sfdefault}
\\begin{document}
\\centering 
"""
    stringwidth=0.08
    

    if subtitle == "":
        baselevel = 1
    else:
        baselevel = 1.5
    stringlist=[[word[0:3],int(word[3:])] for word in string.split('.')]
    stringlist.reverse()
    
    w_current = start
    w_max = start
    for word in stringlist:
        if word[0]=="cap":
            w_current+= 2
            w_max = max(w_max,w_current)
        if word[0]=="cup":
            w_current-= 2
    
    if style=="slices":
        w = w_max+1.5
    else:
        w = w_max
    h = len(stringlist)+baselevel
    dotlength=start
    
    content+="\\psset{{yunit=-1}}\\begin{{pspicture}}(-0.5,-0.5)({0},{1})\n\psset{{linewidth=2.5pt}}\n".format(w-0.5,h+0.5)
    
    content+="\\rput[c]("+str(w/2-0.5)+","+str(0)+"){\\textbf{"+SanitizeTeX(filename)+"}}\n"    
    content+="\\rput[c]("+str(w/2-0.5)+","+str(0.75)+"){"+SanitizeTeX(subtitle)+"}\n"

    def draw_cup(level,index,dotlength):
        string = "\\psbezier({},{})({},{})({},{})({},{})".format(index,level,index,level+0.5,index+1,level+0.5,index+1,level)   
        for i in range(index):
            string += "\\psline({},{})({},{})".format(i,level,i,level+1)
        for i in range(index+2,dotlength):
            string += "\\psbezier({},{})({},{})({},{})({},{})".format(i,level,i,level+0.5,i-2,level+0.5,i-2,level+1)
        return string
    
    def draw_cap(level,index,dotlength):
        string = "\\psbezier({},{})({},{})({},{})({},{})".format(index,level+1,index,level+0.5,index+1,level+0.5,index+1,level+1)
        for i in range(index):
            string += "\\psline({},{})({},{})".format(i,level,i,level+1)
        for i in range(index,dotlength):
            string += "\\psbezier({},{})({},{})({},{})({},{})".format(i,level,i,level+0.5,i+2,level+0.5,i+2,level+1)
        return string
    
    def draw_pos(level,index,dotlength):
        string = "\\psbezier({},{})({},{})({},{})({},{})".format(index,level,index,level+0.5,index+1,level+0.5,index+1,level+1)
        string += "\\psbezier[linecolor=white,linewidth=10pt]({},{})({},{})({},{})({},{})".format(index+1,level,index+1,level+0.5,index,level+0.5,index,level+1)
        string += "\\psbezier({},{})({},{})({},{})({},{})".format(index+1,level,index+1,level+0.5,index,level+0.5,index,level+1)
        for i in range(index):
            string += "\\psline({},{})({},{})".format(i,level,i,level+1)
        for i in range(index+2,dotlength):
            string += "\\psline({},{})({},{})".format(i,level,i,level+1)
        return string

    def draw_neg(level,index,dotlength):
        string = "\\psbezier({},{})({},{})({},{})({},{})".format(index+1,level,index+1,level+0.5,index,level+0.5,index,level+1)
        string += "\\psbezier[linecolor=white,linewidth=10pt]({},{})({},{})({},{})({},{})".format(index,level,index,level+0.5,index+1,level+0.5,index+1,level+1)
        string += "\\psbezier({},{})({},{})({},{})({},{})".format(index,level,index,level+0.5,index+1,level+0.5,index+1,level+1)
        for i in range(index):
            string += "\\psline({},{})({},{})".format(i,level,i,level+1)
        for i in range(index+2,dotlength):
            string += "\\psline({},{})({},{})".format(i,level,i,level+1)
        return string

    # Drawing code
    for level,word in enumerate(stringlist):
        level+=baselevel
        if word[0]=="pos":
            content+=draw_pos(level,word[1],dotlength)
        
        if word[0]=="neg":
            content+=draw_neg(level,word[1],dotlength)
        
        if word[0]=="cup":
            content+=draw_cup(level,word[1],dotlength)
            dotlength-=2
        
        if word[0]=="cap":
            content+=draw_cap(level,word[1],dotlength)
            dotlength+=2
        
        if (style=="slices") and (level>baselevel):
            content+="\\psline[linecolor=lightgray]("+str(w-0.75)+","+str(level)+")(-0.25,"+str(level)+")\n"
        
        if (style=="slices"):
            content+="\\rput[c]("+str(w-1.5)+","+str(level+0.5)+"){\color{gray}"+word[0]+str(word[1])+"}\n"
  
    content+="\\end{pspicture}\n\\end{document}"
    with open("examples/PSTricks/"+filename+"-tangle.tex", "w") as text_file:
        print(content, file=text_file)
        
    run("cd examples/PSTricks && pdflatex -shell-escape "+filename+"-tangle.tex > "+filename+"-tangle.out 2>&1", shell=True)
    run("cd examples/PSTricks && rm "+(" ".join([filename+"-tangle."+string+" " for string in ["log","aux","pdf","out"]])), shell=True)


def drawtangle_old(string,name,style="plain",start=1,subtitle=""):
    """draw the tangle specified by 'string', which is a concatenation of words <type>+<index>, separated by '.' read from right to left, for each elementary tangle slice, read from top to bottom, where:
    <type> is equal to:
        'pos': positive crossing
        'neg': negative crossing
        'cup': cup
        'cap': cap
    <index> is the index at which the crossing, cap or cup sits. 
    The optional parameter 'style' is either 'plain' (simple tangle) or 'slices' (shows slices with labels).
    The optional parameter 'start' is an integer which specifies the number of tangle ends at the top.
    The optional parameter 'subtitle' is a string which is added to the top.
    """
    stringwidth=0.08
    
    if subtitle == "":
        baselevel = 1
    else:
        baselevel = 1.5
    
    # drawing code for elementary slices
    def draw_cup(level,index,dotlength):
        #arcs connecting dots
        ctx.move_to(index,level)
        ctx.curve_to(index,level+0.5,index+1,level+0.5,index+1,level)
               
        for i in range(index):
            ctx.move_to(i,level)
            ctx.line_to(i,level+1)
        
        for i in range(index+2,dotlength):
            ctx.move_to(i,level)
            ctx.curve_to(i,level+0.5,i-2,level+0.5,i-2,level+1)
        
        ctx.set_source_rgb(0,0,0)
        ctx.set_line_width(stringwidth)
        ctx.stroke()
    
    def draw_cap(level,index,dotlength):
        #arcs connecting dots
        ctx.move_to(index,level+1)
        ctx.curve_to(index,level+0.5,index+1,level+0.5,index+1,level+1)
               
        for i in range(index):
            ctx.move_to(i,level)
            ctx.line_to(i,level+1)
        
        for i in range(index,dotlength):
            ctx.move_to(i,level)
            ctx.curve_to(i,level+0.5,i+2,level+0.5,i+2,level+1)
            
        ctx.set_source_rgb(0,0,0)
        ctx.set_line_width(stringwidth)
        ctx.stroke()
    
    def draw_pos(level,index,dotlength):
        #arcs connecting dots
        ctx.move_to(index,level)
        ctx.curve_to(index,level+0.5,index+1,level+0.5,index+1,level+1)
        ctx.set_source_rgb(0,0,0)
        ctx.set_line_width(stringwidth)
        ctx.stroke()
        
        ctx.set_line_cap(cairo.LINE_CAP_BUTT)
        ctx.move_to(index+1,level)
        ctx.curve_to(index+1,level+0.5,index,level+0.5,index,level+1)
        ctx.set_source_rgb(1,1,1)
        ctx.set_line_width(3*stringwidth)
        ctx.stroke()
        
        ctx.set_line_cap(cairo.LINE_CAP_ROUND)
        ctx.move_to(index+1,level)
        ctx.curve_to(index+1,level+0.5,index,level+0.5,index,level+1)
               
        for i in range(index):
            ctx.move_to(i,level)
            ctx.line_to(i,level+1)
        
        for i in range(index+2,dotlength):
            ctx.move_to(i,level)
            ctx.line_to(i,level+1)
        
        ctx.set_source_rgb(0,0,0)
        ctx.set_line_width(stringwidth)
        ctx.stroke()
    
    def draw_neg(level,index,dotlength):
        #arcs connecting dots
        ctx.move_to(index+1,level)
        ctx.curve_to(index+1,level+0.5,index,level+0.5,index,level+1)
        ctx.set_source_rgb(0,0,0)
        ctx.set_line_width(stringwidth)
        ctx.stroke()
        
        ctx.set_line_cap(cairo.LINE_CAP_BUTT)
        ctx.move_to(index,level)
        ctx.curve_to(index,level+0.5,index+1,level+0.5,index+1,level+1)
        ctx.set_source_rgb(1,1,1)
        ctx.set_line_width(3*stringwidth)
        ctx.stroke()
        
        ctx.set_line_cap(cairo.LINE_CAP_ROUND)
        ctx.move_to(index,level)
        ctx.curve_to(index,level+0.5,index+1,level+0.5,index+1,level+1)
               
        for i in range(index):
            ctx.move_to(i,level)
            ctx.line_to(i,level+1)
        
        for i in range(index+2,dotlength):
            ctx.move_to(i,level)
            ctx.line_to(i,level+1)
        
        ctx.set_source_rgb(0,0,0)
        ctx.set_line_width(stringwidth)
        ctx.stroke()
        
    stringlist=[[word[0:3],int(word[3:])] for word in string.split('.')]
    stringlist.reverse()
    scale = 100
    
    w_current = start
    w_max = start
    for word in stringlist:
        if word[0]=="cap":
            w_current+= 2
            w_max = max(w_max,w_current)
        if word[0]=="cup":
            w_current-= 2
    
    if style=="slices":
        w = w_max+1.5
    else:
        w = w_max
    h = len(stringlist)+baselevel
    dotlength=start
    
    surface = cairo.PDFSurface("examples/"+name+'.pdf',w*scale,(h+1)*scale)
    ctx = cairo.Context(surface)
    ctx.set_source_rgb(1, 1, 1) #background colour
    ctx.paint()
    matrix = cairo.Matrix(scale,0,0,scale,0.5*scale,0.5*scale)
    ctx.set_matrix(matrix)
    ctx.set_line_cap(cairo.LINE_CAP_ROUND)
    
    # Drawing code
    for level,word in enumerate(stringlist):
        level+=baselevel
        if word[0]=="pos":
            draw_pos(level,word[1],dotlength)
        
        if word[0]=="neg":
            draw_neg(level,word[1],dotlength)
        
        if word[0]=="cup":
            draw_cup(level,word[1],dotlength)
            dotlength-=2
        
        if word[0]=="cap":
            draw_cap(level,word[1],dotlength)
            dotlength+=2
        
        if (style=="slices") and (level>baselevel):
            ctx.move_to(-0.25,level)
            ctx.line_to(w-0.75,level)
        
            ctx.set_source_rgb(0.8,0.8,0.8)
            ctx.set_line_width(stringwidth)
            ctx.stroke()
        
        if (style=="slices"):
            ctx.set_font_size(0.40)
            ctx.set_source_rgb(0.5,0.5,0.5)
            ctx.select_font_face("Courier",cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_BOLD)
            s = word[0]+str(word[1])
            ctx.move_to(w-2,level+0.5)
            ctx.show_text(s)
    
    s = name
    titlescale = 0.6
    ctx.set_source_rgb(0,0,0)
    ctx.select_font_face("Courier",cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_BOLD)
    tw = w-1
    while tw >= w-1:# adapt font size to width of the picture
        ctx.set_font_size(titlescale)
        tw = ctx.text_extents(s)[2]
        titlescale = titlescale*0.95
    fascent, fdescent, fheight, fxadvance, fyadvance = ctx.font_extents()
    ctx.move_to(w/2-0.5-tw/2.0,fheight/2)
    ctx.show_text(s)  
    
    s = subtitle
    titlescale = 0.4
    ctx.set_source_rgb(0,0,0)
    ctx.select_font_face("Courier",cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_BOLD)
    tw = w-1
    while tw >= w-1:# adapt font size to width of the picture
        ctx.set_font_size(titlescale)
        tw = ctx.text_extents(s)[2]
        titlescale = titlescale*0.95
    fascent, fdescent, fheight, fxadvance, fyadvance = ctx.font_extents()
    ctx.move_to(w/2-0.5-tw/2.0,baselevel-0.75+fheight/2)
    ctx.show_text(s)

################
# Obsolete code:
################

def drawclt(clt,name):
    """Create a pdf file 'name'.pdf in the subfolder 'examples' with a pictographic representation of the CLT 'clt'.
    This function is not used anywhere in the code and kept for debugging purposes only.
    """
    scale = 100
    w = max(clt.top,clt.bot)
    h = max([1]+[abs(i[1]-i[0]) for i in clt.arcs_top()]+\
                 [abs(i[1]-i[0]) for i in clt.arcs_bot()]+\
                 [abs(i[1]-i[0]-clt.top) for i in clt.arcs_mix()])
    
    surface = cairo.PDFSurface("examples/"+name+'.pdf',w*scale,(h+1)*scale)
    ctx = cairo.Context(surface)
    matrix = cairo.Matrix(scale,0,0,scale,0.5*scale,0.5*scale)
    ctx.set_matrix(matrix)
    ctx.set_line_cap(cairo.LINE_CAP_ROUND)
    
    # Drawing code
    draw_tangle_ends(0,0,clt,h,ctx)
    draw_arcs(0,0,clt,h,ctx,(0,0,0))

def drawcob(cob,name):
    """Create a pdf file 'name'.pdf in the subfolder 'examples' with a pictographic representation of the cobordism 'cob'.
    This function is not used anywhere in the code and kept for debugging purposes only.
    """
    clt1=cob.front
    clt2=cob.back
    scale = 100
    w = max(clt1.top,clt1.bot)+2
    h = max([1]+[abs(i[1]-i[0]) for i in clt1.arcs_top()]+\
                 [abs(i[1]-i[0]) for i in clt1.arcs_bot()]+\
                 [abs(i[1]-i[0]-clt1.top) for i in clt1.arcs_mix()]+\
                 [abs(i[1]-i[0]) for i in clt2.arcs_top()]+\
                 [abs(i[1]-i[0]) for i in clt2.arcs_bot()]+\
                 [abs(i[1]-i[0]-clt2.top) for i in clt2.arcs_mix()])
    
    surface = cairo.PDFSurface("examples/"+name+'.pdf',w*scale,(h+1)*scale*len(cob.decos))
    ctx = cairo.Context(surface)
    matrix = cairo.Matrix(scale,0,0,scale,0.5*scale,0.5*scale)
    ctx.set_matrix(matrix)

    ctx.set_line_cap(cairo.LINE_CAP_ROUND)
    
    # Drawing code
    
    for deco_index,deco in enumerate(cob.decos):
        #arcs connecting dots in clt2 blue ("beta")
        draw_arcs(2,(h+1)*deco_index,clt2,h,ctx,(0,0,255))
        #arcs connecting dots in clt1 red ("alpha")
        draw_arcs(2,(h+1)*deco_index,clt1,h,ctx,(255,0,0))
        
        # draw dots on components
        for dot_index,dot in enumerate(deco[1:-1]):
            if dot in [0,1]:
                if dot==1:
                    draw_dot_on_arc(cob.comps[dot_index][0:2],clt1,h,ctx,deco_index)
            else:
                raise Exception('Some components are decorated by more than one dot.')
        
        # draw coefficient + H power
        ctx.set_font_size(0.40)
        ctx.select_font_face("Courier",cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_BOLD)
        ctx.set_source_rgb(0, 0, 0)
        draw_tangle_ends(2,0+(h+1)*deco_index,clt1,h,ctx)
        
        s = '{}.'.format(deco[-1])+'H^{}'.format(deco[0])
        xbearing, ybearing, width, height, dx, dy = ctx.text_extents(s)
        ctx.move_to(0, h/2+(h+1)*deco_index)
        ctx.show_text(s)
    
        ctx.set_source_rgb(0, 0, 0)
        ctx.set_line_width(0.1)
        ctx.stroke()

