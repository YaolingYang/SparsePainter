This example includes 8000 reference individuals from 4 populations with 2091 SNPs (``donor.phase.gz``), and the aim is to paint 500 target individuals (``target.phase.gz``). Remember we have compiled HMPaint in ``HMPaint.exe``, then we can paint with the following command:

``
./HMPaint.exe -reffile donor.phase.gz -targetfile target.phase.gz -popfile popnames.txt -mapfile map.txt -targetname targetname.txt -out HM
``

The output file for this example includes ``HM_painting.txt.gz``, ``HM_aveSNPpainting.txt.gz``, ``HM_aveindpainting.txt.gz`` and ``HM_AAS.txt``.
