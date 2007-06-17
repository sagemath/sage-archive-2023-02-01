"""

We assume a UNIX user named "sageserver".

if .ssh/id_rsa.pub doesn't exist for sageserver user:

 cd
 ssh-keygen -t rsa


# Make these all with same stupid password

write to script

 adduser sage0
 ...


su -c script_name


# set them up for ssh automated logins
 ssh sage1@hostname "mkdir .ssh"
 scp .ssh/id_rsa.pub sage1@hostname:.ssh/authorized_keys2

 ...
# set all passwords of user accounts to be non-login-able.

This should work to run SAGE for all users

(have to automate doing
"Are you sure you want to continue connecting (yes/no)? yes"
)

was@ubuntu:~$ ssh -t sage1@ubuntu sage

# Groups

Change lines in /etc/groups:

  sage1:x:1001:sageserver
  ...
"""
