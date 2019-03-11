from CP2 import *
import sys

def main():
    sirs = SIRS(sys.argv[1],sys.argv[2], sys.argv[3], sys.argv[4])
    init_state = sys.argv[-1]
    if init_state == "Random":
        sirs.initalise_random()
    elif init_state == "Uniform":
        sirs.initalise_uniform()
    else:
        print "Input Error"
        exit()
    Visualise(sirs.update)
main()
