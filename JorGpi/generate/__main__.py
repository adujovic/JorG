import JorGpi.generate.run
from sys import argv

if __name__ == '__main__':
    engine = JorGpi.generate.run.JorGpi(*argv)
    engine.initialize_new_cell()
    engine.possible_configurations()
    exit()
