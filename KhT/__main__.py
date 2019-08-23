import sys

from BNComplexes import *
from CobComplexes import *
from Drawing import *
from CrossingTangle import *

# This is used solely to execute the files in the examples sub directory

def main(data):
    print("----------------")
    print("KhT, version ???")
    print("----------------")
    exec(data)
    print("-------------------------")
    print("KhT executed successfully")
    print("-------------------------")
    
if __name__ == "__main__":
    filename=sys.argv[1]
    try:# assuming KhT calls file without including the '.py'-ending
        with open("examples/"+filename+".py", "r") as text_file:
            data = text_file.read()
    except: # assuming KhT calls file with the '.py'-ending or the file does not have an ending
        with open("examples/"+filename, "r") as text_file:
            data = text_file.read()
    main(data)
