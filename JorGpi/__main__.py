import JorGpi.run
from sys import argv

if __name__ == '__main__':
    engine = JorGpi.run.JorGpi(*argv)
    engine.initialize_new_cell()
    engine.generate_possible_configurations()
    exit()
