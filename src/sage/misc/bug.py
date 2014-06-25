"""
Bug reporting
"""

def bug():
    r"""
    Simple script for putting a bug report in a standard format.

    EXAMPLES:
`       sage: bug()   # not tested
        We will now create a bug report in a standard format.
        Bug in command (one line)? bug
        Expected output (one line)? bugs
        Output you got (one line)? bug report
        Comments or suggestions (one line)? help!
                ******************************
                Please copy and paste the following output into an email (edit
                it further) and send it to sage-support@lists.sourceforge.net.
                ******************************
        From: was
        Machine specs: Linux sha 2.6.16.16 #15 SMP Sun May 14 12:40:14 PDT 2006 i686
        Date: 2006-05-14
         Bug in command (one line): bug
         Exected output: bugs
         Got instead: bug report
         Comments: help!

    AUTHOR: David Joyner (2006-05)
    """
    print("We will now create a bug report in a standard format.")
    import os
    import time
    report = """
        ******************************
        Please copy and paste the following output into an email (edit
        it further) and send it to sage-support@lists.sourceforge.net.
        ******************************
    \n"""
    try:
        lgn = os.getlogin()
    except OSError:
        lgn = os.popen('whoami').read().strip()
    try:
        unm = os.uname()
    except OSError:
        unm = "unknown machine"
    tm = time.strftime('%Y-%m-%d')
    report = report+"From: "+lgn+"\n"
    report = report+"Machine specs: %s\n"%(' '.join(unm))
    report = report+"Date: %s \n"%tm
    report = report + "\n Bug in command (one line): "
    try:
        command = raw_input("Bug in command (one line)? ")
        report = report + command + "\n Exected output: "
        expected = raw_input("Expected output (one line)? ")
        report = report + expected + "\n Got instead: "
        got = raw_input("Output you got (one line)? ")
        report = report + got + "\n Comments: "
        comments = raw_input("Comments or suggestions (one line)? ")
        report = report + comments + "\n"
    except (EOFError, KeyboardInterrupt):
        print("Creation of bug report cancelled.")
        print("Type bug() to create the report again.")
    print(report)

