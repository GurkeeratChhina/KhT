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

import pandas as pd
from graph_tool.all import *
import math
import cairo
from IPython.display import IFrame
from KhT import *
from Tangles import *
from Cobordisms import *
from Complex import *
from tabulate import tabulate

def printdecos(cob,switch="short"):
    if cob.decos==[]:
        return ""
    else:
        if switch == "old long":
            return [cob.comps,cob.decos]
        if switch == "long":
            table=[["H:"]+[deco[0] for deco in cob.decos]]+\
                [[comp]+[deco[i+1] for deco in cob.decos] \
                    for i,comp in enumerate(cob.comps)]+\
                [["coeff:"]+[deco[-1] for deco in cob.decos]]
            tablealt=[["H:",0,1],[[1,2],0,1],["coeff:",3,4]]
            return tabulate(table,tablefmt="plain")
        else:
            return len(cob.decos)
    
def PrettyPrintComplex(Complex,switch="short"): #This is obsolete, use PrettyPrintBNComplex instead
    """Print a complex in human readable form. 
    The second argument is an optional parameter which should be one of the following strings: 
    - 'simple' (default) prints only the length of cobordisms.
    - 'long' prints all cobordism data in a nice table.
    - 'old long' prints all cobordism data, but as a list of lists.
    """
    print("The generators:")
    print(pd.DataFrame({\
        "clt.pairs": [clt.pairs for clt in Complex.elements],\
        "q": [clt.qgr for clt in Complex.elements],\
        "h": [clt.pgr for clt in Complex.elements]
        },columns=["clt.pairs","q","h"]))
    print("The differential: ("+switch+" form)")
    #print(pd.DataFrame([[printdecos(entry,switch)  for entry in row] for row in Complex.morphisms]))
    print(tabulate(pd.DataFrame([[printdecos(entry,switch)  for entry in row] for row in Complex.morphisms]),range(len(Complex.morphisms)),tablefmt="fancy_grid"))
    
def PrintComplexMorphismIntMatrix(Complex): #This is obsolete, use PrettyPrintComplex instead
    for i in Complex.morphisms:
        Row = []
        for j in i:
            Row.append(len(j.decos))
        print(Row)

def PrintComplexMorphismDecoCompMatrix(Complex): #This is obsolete, use PrettyPrintComplex instead
    for i in Complex.morphisms:
        Row = []
        for j in i:
            if j.decos == []:
                Row.append(0)
            else:
                Row.append([j.comps, j.decos])
        print(Row)

def CobordismToDS(Cob): #more pythonic implementation # This is obsolete, done by CobordismToBNAlg
    """ DS is a linear combination of morphisms that are powers of D or powers of S, represented as a list of lists
        an element of DS will be refered to a ds, which is a 3 element list
        the first element is D if the morphism is a power of D, and S if it is a power of S
        the second element is the power of the corresponding morphism
        the third element is the coefficient of the morphism
        Requires Cob to be a cobordism between 4 ended tangles """
    
    if Cob.front.total !=2 or Cob.back.total !=2:
        raise Exception("Cobordism to convert to DS is not between 4-ended tangles")
    DS = [] 
    for elem in Cob.decos:
        if len(elem) == 3: # elem is H^k S
            for ds in DS:
                if ds[0] == "S" and ds[1] == 2*elem[0] +1: #check if saddle is already there
                    ds[2] += elem[2]
                    break
            else:
                DS.append(["S", 2*elem[0]+1, elem[2]]) #otherwise add saddle
        elif len(elem) == 4 and elem[find_first_index(Cob.comps,notcontains_0)+1] == 1: # elem is H^k D
            for ds in DS:
                if ds[0] == "D" and ds[1] == elem[0] +1: #check if dot is already there
                    ds[2] += elem[3]
                    break
            else:
                DS.append(["D", elem[0]+1, elem[3]]) #otherwise add dot
        elif len(elem) == 4 and elem[0] != 0: # elem is H^k times id
            for ds in DS:
                if ds[0] == "S" and ds[1] == 2*elem[0]: #check if saddle part is already there
                    ds[2] += elem[3]
                    break
            else:
                DS.append(["S", 2*elem[0], ((-1)**elem[0])*elem[3]]) #otherwise add saddle part
            for ds in DS:
                if ds[0] == "D" and ds[1] == elem[0]: #check if dot part is already there
                    ds[2] += elem[3]
                    break
            else:
                DS.append(["D", elem[0], elem[3]]) #otherwise add dot part
        else: #elem is id, elem should be [0, 0, 0, n]
            for ds in DS:
                if ds[0] == "S" and ds[1] == 0: #check if id is already there
                    ds[2] += elem[3]
                    break
            else:
                DS.append(["S", 0, elem[3]]) #otherwise add id
    return DS

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

def drawclt(clt,name):
    scale = 100
    w = max(clt.top,clt.bot)
    h = max([1]+[abs(i[1]-i[0]) for i in clt.arcs_top()]+\
                 [abs(i[1]-i[0]) for i in clt.arcs_bot()]+\
                 [abs(i[1]-i[0]-clt.top) for i in clt.arcs_mix()])
    
    surface = cairo.PDFSurface("Output/"+name+'.pdf',w*scale,(h+1)*scale)
    ctx = cairo.Context(surface)
    matrix = cairo.Matrix(scale,0,0,scale,0.5*scale,0.5*scale)
    ctx.set_matrix(matrix)
    ctx.set_line_cap(cairo.LINE_CAP_ROUND)
    
    # Drawing code
    draw_tangle_ends(0,0,clt,h,ctx)
    draw_arcs(0,0,clt,h,ctx,(0,0,0))
    
    return IFrame("Output/"+ name+'.pdf', width='100%', height='300')

def drawcob(cob,name):
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
    
    surface = cairo.PDFSurface("Output/"+name+'.pdf',w*scale,(h+1)*scale*len(cob.decos))
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
    
    return IFrame('Output/' + name+'.pdf', width='100%', height='300')

def DrawFourEndedChainComplex(complex, filename): #TODO: dont call reduce decorations when labeling edges # This is now obsolete, use DrawBNComplex instead
    # For now assume that all CLT are 2-2
    for CLT in complex.elements:
        if CLT.top != 2 or CLT.bot !=2:
            raise Exception("Not a four ended tangle")
    
    # If a CLT is a 2-2 tangle, then the horizontal tangle is [1,0,3,2] and the vertical is [2,3,0,1]
    g = Graph()
    size = len(complex.elements)
    g.add_vertex(size)
    # print("size = " + str(size))
    Vertex_labeling = g.new_vertex_property("string") # "black" is the horizontal CLT and "white" is the vertical CLT
    Position = g.new_vertex_property("vector<float>")
    for i, clt in enumerate(complex.elements):
        if clt.arcs[0] == 1:
            Vertex_labeling[g.vertex(i)] = "black"
        elif clt.arcs[0] == 2:
            Vertex_labeling[g.vertex(i)] = "white"
        else: 
            print(complex.elements[i].arcs)
            print(complex.elements[i].arcs[0])
            raise Exception("CLT to label vertex is not a horizontal or vertical CLT")
    
        #Position[g.vertex(i)] = [50*(2*i+1), 200]
    
    # TODO: omit coefficients of +- 1, leaving only the sign
    SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
    Edge_labeling = g.new_edge_property("string") # construct edge labels with linear combinations of powers of S and D
    for i, j in product(range(size), range(size)): #TODO: rewrite with enumerate
        if complex.morphisms[j][i].ReduceDecorations() != []: #TODO: reduce decorations earlier, not here
            g.add_edge(g.vertex(i), g.vertex(j))
            Edge_labeling[g.edge(i,j)] = ""
            for ds in CobordismToDS(complex.morphisms[j][i]):
                if Edge_labeling[g.edge(i,j)] != "" and ds[2] > 0:
                    Edge_labeling[g.edge(i,j)] += "+"
                if ds[2] != 0:
                    Edge_labeling[g.edge(i,j)] += str(ds[2]) + "·" + ds[0] + str(ds[1]).translate(SUP)

    graph_draw(g, vertex_color = "black", vertex_fill_color = Vertex_labeling, vertex_size = 30,\
                edge_color = "black", edge_pen_width = 6.0, edge_text = Edge_labeling, edge_text_color = "black",\
                edge_text_distance = 10, edge_font_weight = cairo.FONT_WEIGHT_BOLD, edge_font_size = 30, \
                output_size=(200*size, 150*size), edge_marker_size = 20, output="Output/" + filename)
