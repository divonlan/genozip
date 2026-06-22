// ------------------------------------------------------------------
//   text-license.h
//   Copyright (C) 2019-2026 Genozip Limited. Patent Pending.
//   Please see terms and conditions in LICENSE.txt

#include "version.h"

static rom license_text[] = {
    "This program, \"Genozip\", which includes four tools (genozip, genounzip, genocat and genols), source code, object code, executables, documentation and other files, was developed by Divon Lan (\"Developer\") and is copyright (C) 2019-2026 Genozip Limited (\"Licensor\"). All rights reserved. Patent pending.",
    "TERMS AND CONDITIONS FOR USE",
    
    "1. Definitions.",
    "\"License\" shall mean the terms and conditions for use as defined by Sections 1 through 16 of this document.",

    "\"Legal Entity\" shall mean the union of the acting entity and all other entities that control, are controlled by, or are under common control with that entity. For the purposes of this definition, \"control\" means (i) the power, direct or indirect, to cause the direction or management of such entity, whether by contract or otherwise, or (ii) ownership of fifty percent (50%) or more of the outstanding shares, or (iii) beneficial ownership of such entity.",

    "\"You\" (or \"Your\") shall mean Legal Entity (possibly an individual) exercising permissions granted by this License.",
    
    "\"University\" shall mean a Legal Entity that contributes to the scientific record by regularly publishing papers in scientific journals AND that grants academic degrees which are recognized as such by the competent authority in the country in which said Legal Entity is organized.",

    "\"Derivative Works\" shall mean any work that is based on (or derived from) Genozip and for which the editorial revisions, annotations, elaborations, or other modifications represent, as a whole, an original work of authorship. For the purposes of this License, Derivative Works shall not include works that remain separable from Genozip and Derivative Works thereof.",
    
    "\"Your Commercial Data\" shall mean data which You (the Legal Entity exercising permissions granted by this License) obtained with intention of using it in the development process of a product or service and/or was obtained and/or will be used in the context of any kind of service You provide (including also: sequencing, bioinformatics, diagnostic, clinical, research-on-contract, IT, product development, inspection, biobanks and other dataset aggregations for consumption by external users (even if for academic purposes).) for which You get paid, including through fees, grants, salary or government funding. Data derived from Your Commercial Data is also Your Commercial Data.",

    "\"Your Computers\" shall mean computers You own and/or cloud accounts You own at 3rd party cloud providers.",

    "\"Genozip Executables\" shall mean the executable files genozip, genounzip, genocat and genols (with or without an .exe file name suffix).",

    "Other words and terms in this License shall be interpreted as their usual meaning in the context of a software product.",

    "2. Grant of copyright license. Licensor hereby grants to You a limited non-exclusive, non-transferrable, non-sublicensable, revokable copyright license to use Genozip on Your Computers, if you meet the conditions attached to any of the License Types a through g below, for the limited purpose attached to that particular License Type, and subject to the terms and conditions of this License agreement:",
    
    "   a. Standard, Enterprise or Premium License: Using Genozip Executables for any legal purpose, if the license was purchased and paid for, and for the duration that it is in effect. In addition, for Premium License only: Distributing Genozip Executables directly to others (e.g., via email or a private download link). Such distribution does not grant a license to use the software; each recipient must be independently licensed under this Section 2 (such as under a Decompression License) to run it.",

    "   b. Biobank License: Using Genozip Executables to compress data for a public or cross-institutional genomic data repository by certain users named on the license, if the license was purchased and paid for, and for the duration that it is in effect.",
    
    "   c. Research License: Using Genozip Executables for academic research, educational or training purposes provided that You are a University, and limited to a single lab or project within a single University, and the lab or project are not part of a hospital, and the license was purchased and paid for, and for the duration that it is in effect. Use with Your Commercial Data is not permitted.",
    
    "   d. Student License: Using Genozip Executables free of charge, provided You are currently enrolled as a student in an academic degree program at a University, and Your use is strictly limited to Your personal studies. This license automatically expires upon graduation or termination of enrollment, and use with Your Commercial Data is not permitted.",
        
    "   e. Decompression License: Using a subset of Genozip Executables consisting of genounzip, genocat, genols for any legal purpose, on files that were compressed with a valid Genozip license. A Decompression License is free of charge.",
    
    "   f. Evaluation License: Using Genozip Executables for the purpose of evaluating Genozip, free of charge, for a duration of 30 days, if You were not already granted an Evaluation License in the past. The duration of an Evaluation License may be extended by written email approval.",
    
    "   g. Distribution License: For the purpose of hosting and distributing Genozip Executables to others via a public website OR a package or container management system. A Distribution License is free of charge. A Distribution License does not grant a license to use Genozip Executables; each user must be independently licensed under one of the licenses listed in this Section 2.",
    
    "3. Additional Terms and conditions",
    
    "   a. You must fully, truthfully and accurately complete the activation form as prompted by the genozip tool.",
    
    "   b. Using Genozip to compress a file is only permitted if the file is retained in its original form as well or the potential loss of data due to Genozip not being able to uncompress the compressed file would not cause any harm.",
    
    "   c. Any changes to the Genozip's source code and/or creation of Derivative Works (including, but not limited to, forking on github) and/or reverse-engineering of Genozip (except to the extent expressly permitted by applicable mandatory law) and/or using all or part of Genozip's source code (even if modified) in another software package are forbidden, unless prior written permission is obtained from Licensor.",
    
    "   d. Any software source code intentionally submitted for inclusion in Genozip by You to the Licensor or the Developer, including by using a GitHub Pull Request, shall imply complete and irrevocable assignment by You to Licensor of all copyright in the submitted source code. Regarding any such source code You submitted for inclusion in Genozip in the past, You hereby assign all copyright in this submitted source code to Licensor.",
    
    "   e. Reselling Genozip and/or selling a service or a product that includes Genozip or any part of Genozip's code or algorithms (together, \"Genozip Technology\") such that a user of said service or product may directly or indirectly effectuate compression or decompression of data using Genozip Technology - permission for such reselling or selling is not granted in this license, and requires a separate reseller or OEM license. To clarify, merely delivering Genozip-compressed files to others (e.g. your clients or collaborators) IS included in the Standard, Enterprise, Premium, Research and/or Student License and IS NOT subject to this restriction.",

    "   f. Termination and Expiration. This License remains effective until terminated or expired.",
    
    "      (i) Upon Normal Expiration: When Your license expires, concludes, or is no longer valid — including the end of a paid subscription period, the end of a 30-day evaluation, or when a student graduates or ceases active enrollment — Your rights to compress data terminate automatically. However, You retain a free-of-charge Decompression License (as defined in Section 2.e) to uncompress, view, and list files that were compressed with a valid Genozip license, subject to the terms of this License.",

    "      (ii) Termination for Breach: For material breaches as defined in Section 4 termination shall be effective immediately upon such notice. For other breaches, Licensor shall provide written notice specifying the breach, and You shall have fifteen (15) days from receipt of such notice to cure it. If You fail to cure within that period, this License shall terminate automatically. Upon termination for any breach, You must immediately cease all use of Genozip.",

    "   g. Entire Agreement. This License constitutes the entire agreement between You and Licensor relating to Genozip and supersedes all prior or contemporaneous oral or written communications, proposals, representations, and warranties. No amendment to or modification of this License will be binding unless made in a formal written agreement signed by an authorized representative of Licensor.",

    "   h. Authority. The individual accepting this License represents and warrants that they have the authority to bind the Legal Entity identified in the registration. You acknowledge and agree that any acceptance, installation, or use of Genozip by an individual acting on Your behalf shall be binding upon You, regardless of whether that individual had actual authority to act on Your behalf.",

    "4. Remedies for Unauthorized Use. Use of Genozip that is non-compliant with Sections 2, 3a, 3c, or 3e, or use of a Research or Student License to process Your Commercial Data, constitutes a material breach of this License and an infringement of Licensor’s intellectual property rights. In the event of such breach, You agree to: ",

    "   a. Pay Licensor liquidated damages in an amount equal to the greater of: (i) US$100.00 for each file processed during the period of unauthorized use, or (ii) seven (7) times the list price of the applicable commercial license covering the period of unauthorized use, which the parties agree is a reasonable forecast of the commercial value of the processing and not a penalty;",

    "   b. Acknowledge that any data, products, or services developed using the unauthorized software may be deemed infringing, and that Licensor reserves all rights under copyright law to seek injunctive relief to cease the distribution or commercialization of such products or services; and",

    "   c. Reimburse Licensor for all reasonable legal fees, expert fees, and collection costs incurred in enforcing Licensor’s rights under this Section.",

    "5. Data collected. You consent to the following data collection by the Genozip software:",
    
    "   a. At activation time: activation information provided by you and details about your hardware, operating system and IP address.",
    
    "   b. Telemetry: For paid licenses, telemetry collection is opt‑in. For Student and Evaluation licenses, telemetry is always enabled. When a file is compressed, a log record is transmitted containing aggregate statistical information about the compression algorithms' performance and associated metadata. For illustrative examples of the types of data collected, please visit: "WEBSITE_TELEMETRY". The content of that webpage is provided for informational purposes only and does not form part of this License.",
    
    "6. Product Communications. Licensor may send You low-frequency product updates and announcement emails related to Genozip. You may opt out of marketing-specific communications at any time via the provided unsubscribe links, and Your licensing rights are not conditional upon receiving marketing communications.",
    
    "7. Trademarks. This License does not grant permission to use the trade names, trademarks, service marks, or product names of the Licensor, except as required for reasonable and customary use in describing the origin of Genozip.",
    
    "8. Survival. The limitations of liability and ownership rights of Genozip contained herein and Your obligations following termination of this Agreement will survive the termination of this Agreement for any reason.",

    "9. No FDA or other regulatory approvals. The performance characteristics of Genozip have not been established. You acknowledge and agree that (i) Genozip has not been approved, cleared, or licensed by the United States Food and Drug Administration or the Hong Kong Department of Health or any other regulatory entity in any country for any specific intended use, whether research, commercial, diagnostic, or otherwise, and (ii) You must ensure it has any regulatory approvals that are necessary for Your intended uses of Genozip. You will comply with all applicable laws and regulations when using and maintaining Genozip.",
    
    "10. Compliance Verification. Upon reasonable written request by Licensor (no more than once per consecutive twelve (12) month period), You agree to provide a signed statement by an authorized officer or a system-generated usage report certifying that Genozip is being used in full compliance with the terms and licensed scopes of this Agreement. If Licensor reasonably suspects material non-compliance, Licensor reserves the right to have an independent, certified third-party auditor review relevant deployment records under strict confidentiality, during normal business hours, and with at least thirty (30) days' prior written notice.",

    "11. Export Control. You agree to comply with all applicable export control, trade sanctions, and import laws and regulations of the Hong Kong Special Administrative Region and the jurisdiction in which You are incorporated or located. Furthermore, You represent and warrant that You are not located in, organized under the laws of, under the control of, or a national or resident of Iran, Syria, Lebanon, Lybia, or North Korea, nor operating on behalf of any entity domiciled within these territories.",

    "12. Severability. If any provision of this License is held to be unenforceable or invalid, that provision shall be enforced to the maximum extent permissible, and the remaining provisions shall remain in full force and effect.",

    "13. No Waiver. Failure or delay by Licensor in exercising any right under this License shall not constitute a waiver of that right.",

    "14. Governing Law and Jurisdiction. This License shall be governed by and construed in accordance with the laws of the Hong Kong Special Administrative Region. Subject to the final paragraph of this Section, the parties agree that the courts of Hong Kong shall have exclusive jurisdiction to settle any dispute arising under this Agreement. Notwithstanding the foregoing, Licensor (and only Licensor) retains the exclusive right, at its sole discretion, to bring any dispute or claim before an arbitration tribunal administered by the Hong Kong International Arbitration Centre (HKIAC) under its expedited rules, or before any competent court in the jurisdiction where You are incorporated or located.",

    "15. Disclaimer of Warranty. Unless required by applicable law or agreed to in writing, Licensor provides Genozip on an \"AS IS\" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied, including, without limitation, any warranties or conditions of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A PARTICULAR PURPOSE. You are solely responsible for determining the appropriateness of using or redistributing the Genozip and assume any risks associated with Your exercise of permissions under this License.",
    
    "16. LIMITATION OF LIABILITY. TO THE FULLEST EXTENT PERMITTED BY APPLICABLE LAW, IN NO EVENT AND UNDER NO LEGAL THEORY, WHETHER IN TORT (INCLUDING NEGLIGENCE), CONTRACT, STRICT LIABILITY OR OTHER LEGAL OR EQUITABLE THEORY, SHALL LICENSOR OR DEVELOPER BE LIABLE FOR DAMAGES, INCLUDING ANY DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES OF ANY CHARACTER ARISING AS A RESULT OF THIS LICENSE OR OUT OF THE USE OR INABILITY TO USE GENOZIP (INCLUDING BUT NOT LIMITED TO DAMAGES FOR LOSS OF GOODWILL, WORK STOPPAGE, COMPUTER FAILURE OR MALFUNCTION, FILE CORRUPTION, DATA LOSS, OR ANY AND ALL OTHER COMMERCIAL DAMAGES OR LOSSES), EVEN IF LICENSOR OR DEVELOPER HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.  IN NO EVENT WILL LICENSOR'S OR DEVELOPER'S TOTAL LIABILITY TO YOU FOR ALL DAMAGES (OTHER THAN AS MAY BE REQUIRED BY APPLICABLE LAW IN CASES INVOLVING PERSONAL INJURY) EXCEED THE AMOUNT OF $500 USD. THE FOREGOING LIMITATIONS WILL APPLY EVEN IF THE ABOVE STATED REMEDY FAILS OF ITS ESSENTIAL PURPOSE.",
    
    "END OF TERMS AND CONDITIONS", 

    LIC_FIELD_VERSION ": " GENOZIP_CODE_VERSION 
};
 