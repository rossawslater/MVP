from CP2 import *
import sys

def main():
    gol = GoL(int(sys.argv[1]))
    init_state = sys.argv[2]
    if init_state == "Random":
        gol.initalise_random()
    elif init_state == "Uniform":
         gol.insert(gol.bee_hive,int(gol.N/10),int(gol.N/2))
         gol.insert(gol.oscillator,int(gol.N*3/4),int(gol.N/4))
         gol.insert(gol.glider, int(gol.N/2), int(gol.N/2))
    else:
        print "Input Error"
        exit()
    Visualise(gol.update)
main()
