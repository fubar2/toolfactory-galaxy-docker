import argparse
import datetime
import os
import sqlite3



# -g http://localhost:8080 -u $GALAXY_DEFAULT_ADMIN_USER -p $GALAXY_DEFAULT_ADMIN_PASSWORD


def hacktsuser(args):
    dbfile = args.dbfile
    conn = None
    try:
        conn = sqlite3.connect(dbfile)
    except sqlite3.Error as e:
        print(e)

    now = datetime.datetime.now()
    auser = (now, now, "admin@galaxy.org", "fubar", "password", "", "", "", "")
    sql = "UPDATE galaxy_user SET username = 'fubar', password = 'password' WHERE email = 'admin@galaxy.org'"
    #  (create_time,update_time,email,username,password,external,deleted,purged) VALUES(?,?,?,?,?,?,?,?)
    cur = conn.cursor()
    try:
        cur.execute(sql)
        conn.commit()
        print("### executed %s" % sql)
    except:
        print("### failed %s" % sql)
    sql = "INSERT INTO api_keys (user_id,key)  values(1,'fakekey');"
    cur = conn.cursor()
    try:
        cur.execute(sql)
        conn.commit()
        print("### executed %s" % sql)
    except:
        print("### failed %s" % sql)
    sql = "INSERT INTO category (name, description,deleted)  values('Tool Makers','Makes Tools. Is good.',0);"
    cur = conn.cursor()
    try:
        cur.execute(sql)
        conn.commit()
        print("### executed %s" % sql)
    except:
        print("### failed %s" % sql)
    sql = "INSERT INTO category (name, description,deleted)  values('ToolFactory generated tools','ToolFactory products.',0);"
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
    parser.add_argument("-g", "--galaxy", help='URL of target galaxy')
    parser.add_argument("-p", "--password", help='Galaxy admin password')
    parser.add_argument("-e", "--email", help='Galaxy admin email')
    parser.add_argument("-a", "--key", help='Galaxy admin key', default=None)
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

