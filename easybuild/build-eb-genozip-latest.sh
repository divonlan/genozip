#!/bin/bash
genozip_latest_version=`curl -Ls -o /dev/null -w %{url_effective} https://github.com/divonlan/genozip/releases/latest | rev | cut -d- -f1 | rev`
module load EasyBuild
eb --installpath=/apps/skl --try-software-version=$genozip_latest_version genozip.eb
echo "Done building genozip version $genozip_latest_version"
