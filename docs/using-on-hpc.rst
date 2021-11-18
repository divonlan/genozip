.. _using-on-hpc:

Using Genozip on an HPC
=======================

**Building an EasyBuild module**

Please find sample scripts `here <https://github.com/divonlan/genozip/tree/master/easybuild>`_.

**How to register Genozip for a batch job**

Genozip requires registration prior to use, and its a violation of Genozip's license to use Genozip without registration.

`Option 1 (easier)`: register on the login node of the HPC with ``genozip --register``. You need to do this only once.

`Option 2`: If you are able to register on the HPC's login node or your batch script does not have access to your home directory, do the following:

| 1. Register Genozip on another computer with ``genozip --register``. You can skip this step if you have already used (and hence registered) Genozip on this computer. 
| 
| ➔ The license file located at ``~/.genozip_license`` on Linux and Mac and ``%APPDATA%\genozip\.genozip_license`` on Windows.
|
| 2. Copy the license file to the target computer (any directory, any filename). 
|
| 3. Use the ``--licfile`` option to point genozip to the license file, for example:

::

    genozip --licfile mydir/.genozip_license mydata.bam

Still having issues? email support@genozip.com for help.
