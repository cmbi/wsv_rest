# wsv restful service

  This is a restful service that executes whatif wsv commands on request.


## DEPENDENCIES

 * any whatif version that has wsv commands
 * python==3.10.6
 * Flask==2.0.2
 * gevent==21.12.0
 * gunicorn==20.1.0

## INSTALLING

Create a file named wiws/settings.py containing the following:

  WHATIF_EXECUTABLE_PATH = "<path to a script that runs whatif>"
  PDB_DIRECTORY_PATH = "<path to directory containing the rcsb pdb files>"
  PDB_UPLOAD_DIRECTORY_PATH = "<path to directory, for storing uploaded pdb files>"
  SECRET_KEY = "<any random string>"


## RUNNING

To start the server:

  gunicorn -k gevent -b 0.0.0.0:8088 application:app

To access the API, example:

  curl localhost:8088/api/get_hydrogen_bonds_m/101m/
