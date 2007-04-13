Distributed SAGE
Last Updated: 12/14/06
Version Number: 0.1
Contact Info: yqiang@gmail.com

INTRODUCTION

    Distributed SAGE is a program that allows one to easily distribute
    computationally intensive jobs.

OVERVIEW

    There are 3 distinct parts of Distributed SAGE:

    1)  The server

            The server is responsible for job distribution, submission and
            collection.

    2)  The client

            The client is responsible for submitting new jobs to the server
            and collecting the results.

    3)  The worker

            The worker is responsible for the actual processing of the
            jobs.  Workers get jobs from the server, process the jobs and
            then submit the result back to the server.

INSTALLATION

    NEED TO WRITE THIS SECTION


CONFIGURATION

    Before being able to use Distributed SAGE, there are several
    configuration files one must edit.  All configuration files live under
    ~/.dsage.  There are a total of 3 configuration files:

        db.conf (database configuration file)
        server.conf (server configuration file)
        worker.conf (worker configuration file)
        client.conf (client configuration file)

    You will find sample configuration files in distributed_sage/conf. Simply
    copy them to ~./dsage and edit them as neccessary.

    Configuring the server

        1)  Edit server.conf.

        2)  Setting up a public key database
            Authentication of clients is done solely based on public keys.
            To grant a client access to the server, you must first add
            their username and public key to the file you specified in
            server_conf.py under 'pubkey_database'.  The format is simply:

                username:publickey

    Configuring the client

        1)  Edit client.conf.

        2)  Generate a public/private key pair to use for authentication.
            This is something that needs to be done only once.  You want to
            generate a public/private key pair so you can authenticate with
            the server.  To generate a new private/public key, type:

                ssh-keygen -t rsa

            Save the file to a convenient location.  Specify the location of
            your private and public keys in client.conf.  Also, you probably
            want to send your key to the server admin so he can add you to the
            database.

    Configuring the worker
        1)  Edit worker_conf.py
            Leave the id field blank.  It will be automatically generated for
            you the first time you start the worker.


USAGE

    1)  The Distributed SAGE client

        A)  Submitting jobs
            The client submits jobs from the server and awaits results.  To
            start the client type:

                sage console.py

            You will be dropped into a regular python shell.  From there, you
            can proceed as if you were using a normal python shell.
            To instantiate a connection to a server, type:

                P = DSage()

            To submit a job to the server you connected to, type:

                job = P.eval('string to evaluate', 'your job name')

            Optionally, you can also pass the 'preparse=False' keyword
            to eval to stop it from preparsing your string into a SAGE string.

        B)  Checking job results
            To see the result of your computation, simply type:

                print job.result

            It will take a couple of seconds for your result to show up
            even if the computation is trivial, i.e. '2+2' because of the
            overhead involved.  If there is no result yet, you will see:

                'No result yet.'

            Simply check your result later by typing:

                print job.result

        C)  Killing jobs
            To kill a job you submitted, simply do:

                P.kill(job id)

            Note that job id is a integer.  If you do not remember the job
            id, you can simply type:

                print job.result.job.job_id

            This will print out the id of the job you submitted.  Also, the
            DSageInterface() object will keep a list of jobs you have
            submitted in P.jobs. Simply type:

                print P.jobs

            To get a list of jobs.

    2)  The Distributed SAGE worker

        1)  Starting the worker
            To start the worker, simply type:

                    sage worker.py

        2)  Stopping the worker
            To stop the worker, simply hit Ctrl-C.

    3)  The Distributed SAGE server

        A)  Starting the server
            To start the server, simply type:

                sage server.py

        B)  Stopping the server
            To stop the server, simply send it Ctrl-C


VERSION CONTROL
    You can get the latest version of Distributed SAGE by using Mercurial[1].
    The repository is located here:
        http://www.yiqiang.net/dsage/index.cgi/distributed_sage
        http://sage.math.washington.edu/home/yi/Software/index.cgi/distributed_sage
    Use 'hg clone' to make a local copy of the repository, and 'hg pull' to get
    updates.

BUGS
    Send bug reports to yqiang@gmail.com

