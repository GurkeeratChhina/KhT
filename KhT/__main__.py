import sys

import os
from BNComplexes import *
from CobComplexes import *
from Drawing import *
from CrossingTangle import *

# This is used solely to execute the files in the examples subdirectory

def main(data):
    print("----------------")
    print("KhT, version ???")
    print("----------------")
    BNComplexes.filename = filename
    BNComplexes.filepath = filepath
    Drawing.filename = filename
    Drawing.filepath = filepath
    exec(data)
    print("-------------------------")
    print("KhT executed successfully")
    print("-------------------------")
    
if __name__ == "__main__":
    filename=sys.argv[1]
    filepath="examples/"+filename+"/"
    # filename = name of example file called by kht
    # filepath = examples/<filename>/
    
    try:# assuming KhT calls file without including the '.py'-ending
        with open("examples/"+filename+".py", "r") as text_file:
            data = text_file.read()
    except: # assuming KhT calls file with the '.py'-ending or the file does not have an ending
        with open("examples/"+filename, "r") as text_file:
            data = text_file.read()
            filename=filename[:-3]
        
    paths=["","/PSTricks","/BNComplexes","/CobComplexes","Results"]
    for path in paths:
        if not os.path.exists('examples/'+filename+path):
            os.makedirs('examples/'+filename+path)
    
    main(data)
    
    
    
    
