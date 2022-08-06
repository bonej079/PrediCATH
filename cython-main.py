import pyximport;

pyximport.install(language_level=3, inplace=True)

from GoGraph.classes import cliMenu


def main():
    print('Starting the Cli Menu')

    # Fix for Cython path
    import os
    path = os.getcwd() if (os.getcwd().endswith('classes')) else os.path.join(os.getcwd(), 'GoGraph', 'classes')
    print(path)

    cliMenu.runMenu(path)


if __name__ == "__main__":
    main()
