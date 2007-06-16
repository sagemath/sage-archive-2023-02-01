"""

We assume a UNIX user named "sageserver".

if .ssh/id_rsa.pub doesn't exist for sageserver user:

 cd
 ssh-keygen -t rsa


# Make these all with same stupid password
 sudo adduser sage1
 ...

# set them up for ssh automated logins
 ssh sage1@hostname "mkdir .ssh"
 scp .ssh/id_rsa.pub sage1@hostname:.ssh/authorized_keys2

 ...
# set all passwords of user accounts to be random.

This should work to run SAGE for all users

was@ubuntu:~$ ssh -t sage1@ubuntu sage

# Groups

Change lines in /etc/groups:

  sage1:x:1001:sageserver
  ...
"""
