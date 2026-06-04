import sys
from pipeline import Pipeline
from config import Config


if __name__ == "__main__":

    root = sys.argv[1]
    ref = sys.argv[2]

    config = Config()

    Pipeline(root, ref, config).run()