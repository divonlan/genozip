// ------------------------------------------------------------------
//   text-license.h
//   Copyright (C) 2019-2023 Genozip Limited. Patent Pending.
//   Please see terms and conditions in LICENSE.txt

#include "version.h"

static rom license[] = {
    "This program, \"Genozip\", which includes four tools (genozip, genounzip, genocat and genols), source code, object code, executables, documentation and other files, was developed by Divon Lan (\"Developer\") and is copyright (C) 2019-2023 Genozip Limited (\"Licensor\"). All rights reserved. Patent pending.",
    "TERMS AND CONDITIONS FOR USE",
    
    "1. Definitions.",
    "\"License\" shall mean the terms and conditions for use as defined by Sections 1 through 10 of this document.",

    "\"Legal Entity\" shall mean the union of the acting entity and all other entities that control, are controlled by, or are under common control with that entity. For the purposes of this definition, \"control\" means (i) the power, direct or indirect, to cause the direction or management of such entity, whether by contract or otherwise, or (ii) ownership of fifty percent (50%) or more of the outstanding shares, or (iii) beneficial ownership of such entity.",

    "\"You\" (or \"Your\") shall mean Legal Entity (possibly an individual) exercising permissions granted by this License.",

    "\"Derivative Works\" shall mean any work that is based on (or derived from) Genozip and for which the editorial revisions, annotations, elaborations, or other modifications represent, as a whole, an original work of authorship. For the purposes of this License, Derivative Works shall not include works that remain separable from Genozip and Derivative Works thereof.",
    
    "\"Your Commercial Data\" shall mean data which You (the Legal Entity exercising permissions granted by this License) obtained with intention of using it in the development process of a product and/or for provisioning of any kind of service (including also clinical, diagnostic, DNA or RNA sequencing, bioinformatics and cloud services, but excluding education services) for which You get paid. Data derived from Your Commercial Data is also Your Commerical Data.",

    "\"Your Computers\" shall mean computers You own and/or cloud accounts You own at 3rd party cloud providers.",

    "Other words and terms in this License shall be interpreted as their usual meaning in the context of a software product.",

    "2. Grant of copyright license. Licensor hereby grants to You a limited non-exclusive, non-transferrable, non-sublicensable, revokable copyright license to use Genozip on Your Computers, if you meet the conditions attached to any of the License Types a through f below, for the limited purpose attached to that particular License Type, and subject to the terms and conditions of this License agreement:",
    
    "   a. Academic License: For academic research, educational or training purposes provided that You are a recognized academic research institution which is not a hospital, or a registered student at such an institution, but excluding use with Your Commercial Data. An Academic License is free of charge.",
    
    "   b. Academic License: For another non-commercial purpose, if it has been pre-approved by Licensor in writing. Email "EMAIL_REGISTER" to seek such an approval.",
    
    "   c. Standard, Enterprise or Premium License: For any legal purpose, if the license was purchased and paid for, and for the duration that it is in effect.",
    
    "   d. Decompression License: For any legal purpose, if use is limited to the three tools: genounzip, genocat, genols. A Decompression License is free of charge.",
    
    "   e. Evaluation License: For the purpose of evaluating Genozip, free of charge, for a duration of 30 days, if You were not already granted an Evaluation License in the past. The duration of an Evaluation License may be extended by written email approval.",
    
    "   f. Distribution License: For the purpose of distributing Genozip to others via a platform that is free of charge, including (but not limited to) an Internet website, a package or container management system, or a module on an institutional HPC. A Distribution License is free of charge. Each end user must independently register to Genozip and be granted a Standard, Enterprise, Premium, Academic, Decompression or Evaluation License.",
    
    "3. Additional Terms and conditions",
    
    "   a. You must fully, truthfully and accurately complete the registration, either by completing the registration as prompted by the genozip tool or by receiving registration confirmation after registering by emailing "EMAIL_REGISTER". The only exception is if you are registering for the purpose of peer-reviewing a manuscript submitted to a journal by the Genozip author - in this case only, you may enter arbitrary information to preserve anonymity.",
    
    "   b. Using Genozip to compress a file is only permitted if the file is retained in its original form as well or the potential loss of data due to Genozip not being able to uncompress the compressed file would not cause any harm.",
    
    "   c. Any changes to the Genozip's source code and/or creation of Derivative Works and/or reverse-engineering of Genozip and/or using all or part of Genozip's source code (even if modified) in another software package are forbidden, unless prior written permission is obtained from Licensor.",
    
    "   d. Any software source code intentionally submitted for inclusion in Genozip by You to the Licensor or the Developer, including by using a Github Pull Request, shall imply complete and irrevocable assignment by You to Licensor of all copyright in the submitted source code. Regarding any such source code You submitted for inclusion in Genozip in the past, You hereby assign all copyright in this submitted source code to Licensor.",
    
    "   e. Reselling Genozip and/or selling a service or a product that includes Genozip or any part of Genozip's code or algorithms (together, \"Genozip Technology\") such that a user of said service or product may directly or indirectly effectuate compression or decompression of data using Genozip Technology - permission for such reselling or selling is not granted in this license, and requires a separate reseller or OEM license. To clarify, merely delivering Genozip-compressed files to others (e.g. your clients or collaborators) IS included in the Standard, Enterprise, Premium and/or Academic License and IS NOT subject to this restriction.",

    "4. Severely Unauthorized Use of Genozip. Use which is either: a. with no license (Standard License, Enterprise License, Premium License, Acedemic License, Evaluation License, Decompression License or Distribution License) granted according to section 2; or b. non-compliant with any of sections 3a, 3c, 3e; or both - shall be considered Severely Unauthorized Use of Genozip. In this case, You agree to a. pay Licensor US$100.00 for each file You compressed with Genozip b. that Licensor shall be eligible for 20% ownership of any revenue generated and intellectual property created that involved the Severely Unauthorized Use of Genozip. c. reimburse licensor for all legal and/or collection costs related to Your Severely Unauthorized Use of Genozip.",

    "5. Data collected. You consent to the following data collection:",
    
    "   a. At registration time: registration information provided by you and details about your hardware, operating system and IP address as displayed at end of the registration process.",
    
    "   b. When a file is compressed: a log record containing aggregate statistical information about the performance of the compression algorithm and associated metadata. This logging occurs when using an Academic or Evaluation license. When using a Standard, Enterprise or Premium (i.e. paid) License, YOU WILL BE ASKED TO CHOOSE WHETHER OR NOT YOU ALLOW THIS LOGGING. Details can be found here: " WEBSITE_LOGS". ",
 
    "6. Trademarks. This License does not grant permission to use the trade names, trademarks, service marks, or product names of the Licensor, except as required for reasonable and customary use in describing the origin of the Genozip",
    
    "7. Survival. The limitations of liability and ownership rights of Genozip contained herein and Licensee’s obligations following termination of this Agreement will survive the termination of this Agreement for any reason.",

    "8. No FDA or other regulatory approvals. The performance characteristics of Genozip have not been established. Licensee acknowledges and agrees that (i) Genozip has not been approved, cleared, or licensed by the United States Food and Drug Administration or the Hong Kong Department of Health or any other regulatory entity in any country for any specific intended use, whether research, commercial, diagnostic, or otherwise, and (ii) Licensee must ensure it has any regulatory approvals that are necessary for Licensee’s intended uses of Genozip. Licensee will comply with all applicable laws and regulations when using and maintaining Genozip.",
    
    "9. Disclaimer of Warranty. Unless required by applicable law or agreed to in writing, Licensor provides Genozip on an \"AS IS\" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied, including, without limitation, any warranties or conditions of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A PARTICULAR PURPOSE. You are solely responsible for determining the appropriateness of using or redistributing the Genozip and assume any risks associated with Your exercise of permissions under this License.",
    
    "10. LIMITATION OF LIABILITY. TO THE FULLEST EXTENT PERMITTED BY APPLICABLE LAW, IN NO EVENT AND UNDER NO LEGAL THEORY, WHETHER IN TORT (INCLUDING NEGLIGENCE), CONTRACT, STRICT LIABILITY OR OTHER LEGAL OR EQUITABLE THEORY, SHALL LICENSOR OR DEVELOPER BE LIABLE FOR DAMAGES, INCLUDING ANY DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES OF ANY CHARACTER ARISING AS A RESULT OF THIS LICENSE OR OUT OF THE USE OR INABILITY TO USE GENOZIP (INCLUDING BUT NOT LIMITED TO DAMAGES FOR LOSS OF GOODWILL, WORK STOPPAGE, COMPUTER FAILURE OR MALFUNCTION, FILE CORRUPTION, DATA LOSS, OR ANY AND ALL OTHER COMMERCIAL DAMAGES OR LOSSES), EVEN IF LICENSOR OR DEVELOPER HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.  IN NO EVENT WILL LICENSOR'S OR DEVELOPER'S TOTAL LIABILITY TO LICENSEE FOR ALL DAMAGES (OTHER THAN AS MAY BE REQUIRED BY APPLICABLE LAW IN CASES INVOLVING PERSONAL INJURY) EXCEED THE AMOUNT OF $500 USD. THE FOREGOING LIMITATIONS WILL APPLY EVEN IF THE ABOVE STATED REMEDY FAILS OF ITS ESSENTIAL PURPOSE.",
    
    "END OF TERMS AND CONDITIONS",

    LIC_FIELD_VERSION ": " GENOZIP_CODE_VERSION 
};
 