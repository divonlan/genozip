.. _using-on-hpc:

I can't register Genozip's license
==================================

Genozip requires registration prior to use, and its a violation of Genozip's non-commercial license to use Genozip without registration.

However, in some environments, interactive registration is not possible. For example, when submitting a batch job to an HPC.

In these environments, do the following:

| 1. Register Genozip on another computer. This creates a license file located at ``~/.genozip_license`` on Linux and Mac and ``%APPDATA%/genozip/.genozip_license`` on Windows.
|
| 2. Copy the license file to the target computer (any directory, any filename). 
|
| 3. Use the --licfile option to point genozip to the license file, for example:

::

    genozip --licfile mydir/.genozip_license mydata.bam
