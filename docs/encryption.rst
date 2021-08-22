.. _encryption:

Encryption
==========

Many organizations and research projects need to comply with security and privacy requirements in regards to their genomic data.

These days, where news about system security breaches has become routine, it is simply not acceptable anymore to rely solely on the system defenses against intruders, and it is critical to protect the data itself too.

Genozip offers an easy way to encrypt the data while compressing it. The data is encrypted using the `Advanced Encryption Standard (AES) <https://en.wikipedia.org/wiki/Advanced_Encryption_Standard>`_ established by the U.S. National Institute of Standards and Technology (NIST). This is the same encryption method routinely used to protect financial transactions and other sensitive data. Genozip uses the strongest version of AES - 256 bits.

Encryption and decryption are done by simply providing a password during compression and decompression, for example:

``genozip --password mysecret000@ myfile.bam``

``genounzip --password mysecret000@ myfile.bam.genozip``

Behind the scenes, Genozip uses the password to generate an AES encryption key, that is then used to encrypt the data during compression or decrypt it during decompression. The password, and the key derived from it, are not stored anywhere, and Genozip has no way to recover them if they are lost. This means that even if someone successfully breaks into your computer, they still cannot read your data.

