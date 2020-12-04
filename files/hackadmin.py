import argparse
import datetime
import os
import sqlite3



# -g http://localhost:8080 -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD

dome=[ "INSERT INTO galaxy_user VALUES(1,'2020-12-04 02:56:35.066506','2020-12-04 02:56:35.066519','admin@galaxy.org','fubar','fakekey',0,0,0,0)",
 "UPDATE galaxy_user SET username = 'fubar', password = '5baa61e4c9b93f3f0682250b6cf8331b7ee68fd8' WHERE email = 'admin@galaxy.org'",
  "INSERT INTO role VALUES(2,'2020-12-04 02:56:35.083902','2020-12-04 02:56:35.083909','admin@galaxy.org','Private Role for admin@galaxy.org','private',0)",
  "INSERT INTO user_role_association VALUES(1,1,2,'2020-12-04 02:56:35.097627','2020-12-04 02:56:35.097634')",
 "INSERT INTO api_keys (user_id,key)  values(1,'fakekey');",
 "INSERT INTO category (name, description,deleted)  values('Tool Makers','Makes Tools. Is good.',0);",
 "INSERT INTO category (name, description,deleted)  values('ToolFactory generated tools','ToolFactory products.',0);",
]

def hacktsuser(args):
    dbfile = args.dbfile
    conn = None
    try:
        conn = sqlite3.connect(dbfile)
    except sqlite3.Error as e:
        print(e)

    now = datetime.datetime.now()
    for sql in dome:
        cur = conn.cursor()
        try:
            cur.execute(sql)
            conn.commit()
            print("### executed %s" % sql)
        except:
            print("### failed %s" % sql)
    conn.close()

def _parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dbfile", default="/galaxy-central/database/community.sqlite",help='Path to toolshed sqlite db')
    return parser

def main():
    """
    load a folder of histories or a single gz
    """
    args = _parser().parse_args()
    hacktsuser(args)

if __name__ == "__main__":
    main()

