import argparse
import os


from bioblend import galaxy


def _parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--galaxy", help='URL of target galaxy')
    parser.add_argument("-p", "--password", help='Galaxy admin password')
    parser.add_argument("-e", "--email", help='Galaxy admin email')
    parser.add_argument("-a", "--key", help='Galaxy admin key', default=None)
    parser.add_argument("-t", "--tool_id", help='Install conda deps for this tool id')
    return parser

def main():
    """
    install_dependencies(self, tool_id)
    """
    args = _parser().parse_args()
    if args.key:
        gi = galaxy.GalaxyInstance(url=args.galaxy, key=args.key)
    else:
        gi = galaxy.GalaxyInstance(url=args.galaxy, email=args.email, password=args.password)
    x = gi.tools.install_dependencies(tool_id = args.tool_id.strip())
    print(f"Called install_dependencies on {args.tool_id} - got {x}")


if __name__ == "__main__":
    main()

