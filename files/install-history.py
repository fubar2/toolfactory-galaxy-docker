import argparse
import os


from bioblend import galaxy


def _parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--galaxy", help='URL of target galaxy')
    parser.add_argument("-p", "--password", help='Galaxy admin password')
    parser.add_argument("-e", "--email", help='Galaxy admin email')
    parser.add_argument("-a", "--key", help='Galaxy admin key', default=None)
    parser.add_argument("-i", "--history_path", help='Path to history gz files to be loaded')
    return parser

def main():
    """
    load a folder of histories or a single gz
    """
    args = _parser().parse_args()
    if args.key:
        gi = galaxy.GalaxyInstance(url=args.galaxy, key=args.key)
    else:
        gi = galaxy.GalaxyInstance(url=args.galaxy, email=args.email, password=args.password)
    hdir = args.history_path
    # h = gi.histories.get_most_recently_used_history()
    if os.path.isdir(hdir):
        for fp in os.listdir(hdir):
            hp = os.path.join(hdir,fp)
            if os.path.isfile(hp):
                x = gi.histories.import_history(file_path=hp, url=None)
                print('installed ',hp,'res=',x)
    else:
        x = gi.histories.import_history(file_path=hdir, url=None)
        print('installed',hdir,'res=',x)


if __name__ == "__main__":
    main()

